!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for momentum matrix elements, see doi:
!! 10.1103/PhysRevB.63.201101 10.1103/PhysRevB.72.125105 10.1103/PhysRevB.98.115115
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

  !> Calculate momentum matrix elements
  subroutine momentumMatrix(env, parallelKS, eigvals, filling, eigVecsCplx, ham, over,&
      & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, kPoint, kWeight,&
      & cellVec, iCellVec)

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

    !> Sparse Overlap
    real(dp), intent(in) :: over(:)

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

    integer :: iK, nTransition, iKS, iS, iCart, jCart, nOrbs, nSpin, nKpts, nIndepHam, ii, jj
    integer :: iTrans, jTrans, kTrans
    real(dp) :: maxFill
    complex(dp), allocatable :: work(:,:), work2(:,:), dipole(:,:), chi(:,:,:)
    integer, allocatable :: startingState(:,:), endingState(:,:), indx(:)
    real(dp), allocatable :: transitionEgy(:)
    type(TFileDescr) :: fd1

  #:if not WITH_SCALAPACK

    nOrbs = size(eigVals,dim=1)
    nSpin = size(eigVals,dim=3)
    nKpts = size(eigVals,dim=2)

    nIndepHam = nIndependentHam(nSpin)
    maxFill = maximumFillings(nSpin)
    call filledOrEmpty(startingState, endingState, nIndepHam, nKpts, filling, nOrbs, maxFill)

    ! Count total number of transitions at all k-points
    nTransition = sum(startingState*(nOrbs-endingState+1))
    allocate(transitionEgy(nTransition), source = 0.0_dp)
    allocate(chi(3,3,nTransition), source = cmplx(0,0,dp))
    allocate(indx(nTransition))

    ! largest number of transitions for any k-point
    nTransition = maxval(startingState)*(nOrbs-minval(endingState)+1)
    allocate(dipole(nTransition, 3), source = cmplx(0,0,dp))

    allocate(work(nOrbs, nOrbs))
    allocate(work2(nOrbs, nOrbs))

    ! evaluate energy differences for each transition in the system
    nTransition = 0
    do iK = 1, nKpts
      do iS = 1, nIndepHam
        do ii = 1, startingState(iS, iK)
          do jj = endingState(iS, iK), nOrbs
            nTransition = nTransition + 1
            transitionEgy(nTransition) = eigvals(jj, iK, iS) - eigvals(ii, iK, iS)
          end do
        end do
      end do
    end do

    call index_heap_sort(indx, transitionEgy)

    iTrans = 0
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      do iCart = 1, 3

        ! First term of (7) from 10.1103/PhysRevB.98.115115 since basis is non-orthogonal
        work2(:,:) = cmplx(0,0,dp)
        call unpackHSdk(work2, over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)
        call hemm(work, 'l', work2, eigVecsCplx(:,:,iKS))
        do ii = 1, nOrbs
          work(:,ii) = -eigvals(ii, iK, iS) * work(:,ii)
        end do
        work2(:,:) = cmplx(0,0,dp)
        call unpackHSdk(work2, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)
        call hemm(work, 'l', work2, eigVecsCplx(:,:,iKS), beta=(1.0_dp,0.0_dp))

        jTrans = 0
        do ii = 1, startingState(iS, iK)
          do jj = endingState(iS, iK), nOrbs
            jTrans = jTrans + 1
            dipole(jTrans, iCart) = dot_product(eigVecsCplx(:,jj,iKS), work(:,ii))
          end do
        end do

      end do

      jTrans = 0
      do ii = 1, startingState(iS, iK)
        do jj = endingState(iS, iK), nOrbs
          jTrans = jTrans + 1
          kTrans = iTrans + jTrans
          do iCart = 1, 3
            do jCart = 1, 3
              chi(jCart, iCart, kTrans) = kweight(iK) * dipole(jTrans, iCart)&
                  & * conjg(dipole(jTrans, jCart)) / (pi * transitionEgy(kTrans)**2)
            end do
          end do
        end do
      end do

      iTrans = iTrans + startingState(iS,iK)*(nOrbs-endingState(iS,iK)+1)

    end do

    call openFile(fd1, "dielectric.dat", mode="w")

    do ii = 1, nTransition
      jj = indx(ii)
      write(fd1%unit, "(F10.6,E16.8)")transitionEgy(jj) * Hartree__eV, real(chi(1,1,jj))
    end do

    call closeFile(fd1)

  #:endif

  end subroutine momentumMatrix


!  !> Indexing for transition based on start and end transitions
!  pure function indx2transition(iStart, iEnd, iKS, localKS, startingStates, endingStates, nOrbs)&
!      & result(iTrans)
!
!    integer, intent(in) :: iStart, iEnd, iKS, nOrbs
!    integer, intent(in) :: localKS(:,:), startingStates(:,:), endingStates(:,:)
!
!    integer :: iTrans, jKS, iS, iK
!
!    iTrans = 0
!    do jKS = 1, iKS -1
!      iK = localKS(1, jKS)
!      iS = localKS(2, jKS)
!      iTrans = iTrans + startingStates(iS,iK)*(nOrbs-endingStates(iS,iK)+1)
!    end do
!    iK = localKS(1, iKS)
!    iS = localKS(2, iKS)
!    iTrans = iTrans + (iStart-1) * (nOrbs-endingStates(iS, iK)+1) + iEnd - endingStates(iS, iK) + 1
!
!  end function indx2transition

end module dftbp_derivs_matrixelements
