!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for momentum matrix elements, see doi:
!! 10.1103/PhysRevB.63.201101  10.1103/PhysRevB.98.115115
module dftbp_derivs_matrixelements
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : Hartree__eV, imag, pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_derivs_fillings, only : filledOrEmpty, nIndependentHam, maximumFillings
  use dftbp_dftb_periodic, only : TNeighbourList, TSymNeighbourList
  use dftbp_dftb_nonscc, only : buildS
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
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
  public :: imDielectric, momentumMatrix, approxAtomDipole

contains


  subroutine imDielectric(cellVol, filling, kWeight, pMatrixElements)

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> Occupations of orbitals
    real(dp), intent(in) :: filling(:,:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    complex(dp), intent(in) :: pMatrixElements(:,:)


  end subroutine imDielectric


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
    complex(dp), allocatable :: work(:,:), work2(:,:), pMomentum(:,:)
    integer, allocatable :: nFilled(:,:), nEmpty(:,:), indx(:)
    real(dp), allocatable :: transitionEgy(:), chi(:,:,:)
    type(TFileDescr) :: fd1

  #:if not WITH_SCALAPACK

    nOrbs = size(eigVals,dim=1)
    nSpin = size(eigVals,dim=3)
    nKpts = size(eigVals,dim=2)

    nIndepHam = nIndependentHam(nSpin)
    maxFill = maximumFillings(nSpin)
    call filledOrEmpty(nFilled, nEmpty, nIndepHam, nKpts, filling, nOrbs, maxFill)

    ! Count total number of transitions at all k-points
    nTransition = sum(nFilled*(nOrbs-nEmpty+1))
    allocate(transitionEgy(nTransition), source = 0.0_dp)
    allocate(chi(3,3,nTransition), source = 0.0_dp)
    allocate(indx(nTransition))

    ! largest number of transitions for any k-point
    !nTransition = maxval(nFilled*(nOrbs-nEmpty))+1
    nTransition = maxval(nFilled)*(nOrbs-minval(nEmpty)+1)
    allocate(pMomentum(nTransition, 3), source = cmplx(0,0,dp))

    allocate(work(nOrbs, nOrbs))
    allocate(work2(nOrbs, nOrbs))

    ! evaluate energy differences for each transition in the system
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
        call hemm(work, 'l', work2, eigVecsCplx(:,:,iKS), beta=cmplx(1,0,dp))

        jTrans = 0
        do ii = 1, nFilled(iS, iK)
          do jj = nEmpty(iS, iK), nOrbs
            jTrans = jTrans + 1
            pMomentum(jTrans, iCart) = dot_product(eigVecsCplx(:,jj,iKS), work(:,ii))
          end do
        end do

      end do

      jTrans = 0
      do ii = 1, nFilled(iS, iK)
        do jj = nEmpty(iS, iK), nOrbs
          jTrans = jTrans + 1
          kTrans = iTrans + jTrans
          do iCart = 1, 3
            do jCart = 1, 3
              chi(jCart, iCart, kTrans) = kweight(iK) * real(pMomentum(jTrans, iCart)&
                  & * conjg(pMomentum(jTrans, jCart)), dp)
            end do
          end do
        end do
      end do

      iTrans = iTrans + nFilled(iS,iK)*(nOrbs-nEmpty(iS,iK)+1)

    end do

    call openFile(fd1, "dielectric.dat", mode="w")

    do ii = 1, nTransition
      jj = indx(ii)
      write(fd1%unit, "(F10.6,E16.8)")transitionEgy(jj) * Hartree__eV, chi(1,1,jj)
    end do

    call closeFile(fd1)

  #:endif

  end subroutine momentumMatrix


  !> Evaluate onsite dipole matrix elements using the approximation in Sandu PRB, doi:
  !> 10.1103/PhysRevB.72.125105
  subroutine approxAtomDipole(over, nNeighbour, iNeighbour, iSparseStart, img2CentCell, orb,&
      & species, coord)

    !!> On-site dipole matrix elements to second order in overlap
    !real(dp), intent(out) :: dab(:,:,:,:)

    !> Overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours for overlap for each of the central cell atoms
    integer, intent(in) :: nNeighbour(:)

    !> Neighbour index
    integer, intent(in) :: iNeighbour(0:,:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species-list of atoms
    integer, intent(in) :: species(:)

    !> List of all atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    real(dp) :: tmp(orb%mOrb,orb%mOrb), tmpA(orb%mOrb,orb%mOrb)
    integer :: nAtom0, iAt1, iAt2, iAt2f, iNeigh
    integer :: iSp1, iSp2, nOrb1, nOrb2, iOrig, iCart

    real(dp) :: dab(orb%mOrb,orb%mOrb,size(nNeighbour),3)

    nAtom0 = size(nNeighbour)
    @:ASSERT(all(shape(dab) == [orb%mOrb, orb%mOrb, nAtom0, 3]))
    dab(:,:,:,:) = 0.0_dp

    do iAt1 = 1, nAtom0
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 0, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        nOrb2 = orb%nOrbSpecies(iSp2)
        iOrig = iSparseStart(iNeigh, iAt1) + 1

        tmp(:nOrb2, :nOrb1) = reshape(over(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2, nOrb1])

        tmpA(:nOrb1,:nOrb1) = matmul(transpose(tmp(:nOrb2, :nOrb1)), tmp(:nOrb2, :nOrb1))
        do iCart = 1, 3
          dab(:nOrb1, :nOrb1, iAt1, iCart) = dab(:nOrb1, :nOrb1, iAt1, iCart)&
              & +(coord(iCart, iAt2) - coord(iCart, iAt1)) * tmpA(:nOrb1,:nOrb1)
        end do

        tmpA(:nOrb2,:nOrb2) = matmul(tmp(:nOrb2, :nOrb1), transpose(tmp(:nOrb2, :nOrb1)))
        do iCart = 1, 3
          dab(:nOrb2, :nOrb2, iAt2f, iCart) = dab(:nOrb2, :nOrb2, iAt2f, iCart)&
              & +(coord(iCart, iAt1) - coord(iCart, iAt2)) * tmpA(:nOrb2,:nOrb2)
        end do

      end do
    end do

    dab(:,:,:,:) = 0.25_dp * dab

    do iAt1 = 1, nAtom0
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      write(*,*)iAt1
      do iCart = 1, 3
        do iOrig = 1, nOrb1
          write(*,*)dab(:nOrb1, iOrig, iAt1, iCart)
        end do
        write(*,*)
      end do
    end do

  end subroutine approxAtomDipole

end module dftbp_derivs_matrixelements
