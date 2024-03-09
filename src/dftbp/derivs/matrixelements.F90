!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for momentum matrix elements, see https://doi.org/10.1103/PhysRevB.98.115115
module dftbp_derivs_matrixelements
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : Hartree__eV, imag, pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_derivs_fillings, only : filledOrEmpty, nIndependentHam, maximumFillings
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_parallelks, only : TParallelKS
#:if WITH_SCALAPACK
  use dftbp_dftb_sparse2dense, only : unpackHSdk
  use dftbp_extlibs_scalapackfx, only : CSRC_, DLEN_, MB_, NB_, RSRC_, pblasfx_phemm,&
      & scalafx_getdescriptor
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
#:else
  use dftbp_dftb_sparse2dense, only : unpackHSdk
  use dftbp_math_blasroutines, only : hemm
#:endif

  implicit none

  private
  public :: momentumMatrix

contains

  !> Calculate derivatives with respect to k in periodic structures
  subroutine momentumMatrix(env, parallelKS, eigvals, filling, eigVecsCplx, ham, neighbourList,&
      & nNeighbourSK, denseDesc, iSparseStart, img2CentCell, kPoint, kWeight, cellVec, iCellVec)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Occupations of orbitals
    real(dp), intent(in) :: filling(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    integer :: iK, nTransition, iKS, iS, iCart, nOrbs, nSpin, nKpts, nIndepHam, ii, jj
    real(dp) :: maxFill
    complex(dp), allocatable :: dHamSqr(:,:), work(:,:), dipole(:,:), chi(:,:)
    integer, allocatable :: nFilled(:,:), nEmpty(:,:), indx(:)
    real(dp), allocatable :: transitionEgy(:)
    type(TFileDescr) :: fd1

    nOrbs = size(eigVals,dim=1)
    nSpin = size(eigVals,dim=3)
    nKpts = size(eigVals,dim=2)

    nIndepHam = nIndependentHam(nSpin)
    maxFill = maximumFillings(nSpin)
    call filledOrEmpty(nFilled, nEmpty, nIndepHam, nKpts, filling, nOrbs, maxFill)

    nTransition = 0
    do iK = 1, nKpts
      do iS = 1, nIndepHam
        do ii = 1, nFilled(iS, iK)
          do jj = nEmpty(iS, iK), nOrbs
            nTransition = nTransition + 1
          end do
        end do
      end do
    end do

    allocate(transitionEgy(nTransition), source = 0.0_dp)
    allocate(dipole(nTransition, 3), source = cmplx(0,0,dp))
    allocate(chi(3,nTransition), source = cmplx(0,0,dp))
    allocate(indx(nTransition))

 #:if WITH_SCALAPACK

    do iCart = 1, 3
      dHamSqr(:,:) = 0.0_dp
      call unpackHSdk(env%blacs, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, dHamSqr,&
          & iCart)

    end do

  #:else

    allocate(dHamSqr(nOrbs, nOrbs))
    allocate(work(nOrbs, nOrbs))

    nTransition = 0
    do iK = 1, nKpts
      do iS = 1, nIndepHam
        do ii = 1, nFilled(iS, iK)
          do jj = nEmpty(iS, iK), nOrbs
            nTransition = nTransition + 1
            transitionEgy(nTransition) = eigvals(jj, iK, iS) - eigvals(ii, iK, iS)
          end do
        end do
      end do
    end do

    call index_heap_sort(indx, transitionEgy)

    do iCart = 1, 3

      nTransition = 0
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)

        dHamSqr(:,:) = cmplx(0,0,dp)
        call unpackHSdk(dHamSqr, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)

        call hemm(work, 'l', dHamSqr, eigVecsCplx(:,:,iKS))

        do ii = 1, nFilled(iS, iK)
          do jj = nEmpty(iS, iK), nOrbs
            nTransition = nTransition + 1
            dipole(nTransition, iCart) = dot_product(eigVecsCplx(:,jj,iKS), work(:,ii))
            chi(iCart, nTransition) = pi * kWeight(iK) * &
                & (filling(ii,iK,iS) - filling(jj,iK,iS)) *&
                & (dipole(nTransition, iCart))**2 / transitionEgy(nTransition)**2
          end do
        end do
      end do

    end do

    call openFile(fd1, "dielectric.dat", mode="w")

    do ii = 1, nTransition
      jj = indx(ii)
      write(fd1%unit,*)transitionEgy(jj) * Hartree__eV, abs(sum(chi(:,jj)))/3.0_dp
    end do

    call closeFile(fd1)

  #:endif

  end subroutine momentumMatrix

end module dftbp_derivs_matrixelements
