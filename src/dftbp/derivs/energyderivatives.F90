!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module containing routines for energy derivatives, based on expressions from Gonze and Vigneron,
!> Phys Rev B 39 13120 (1989).
module dftbp_derivs_energyderivatives
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_shift, only : addShift, totalShift
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  implicit none

  private
  public :: secondEDerivative

  interface secondEDerivative
    module procedure secondEDerivative_real
  end interface secondEDerivative

contains

  !> Energy second derivatives
  subroutine secondEDerivative_real(env, d2E, dPsiReal1, dPsiReal2, dPotential1, dPotential2,&
      & d2Potential, eigvals, filling, eigVecsReal, nFilled, over, ham, nNeighbourSK,&
      & neighbourList, species, orb, denseDesc, iSparseStart, nAtom, img2CentCell)

    !> Computational environment settings
    type(TEnvironment), intent(inout) :: env

    !> Second derivative of the energy
    real(dp), intent(out) :: d2E

    !> Derivatives of wavefunctions with respect to perturbations (nOrb, nOrb)
    real(dp), intent(in) :: dPsiReal1(:,:)

    !> Derivatives of wavefunctions with respect to perturbations (nOrb, nOrb)
    real(dp), intent(in) :: dPsiReal2(:,:)

    !> 1st derivatives of external potential (mOrb, mOrb, nAtom, nSpin)
    real(dp), intent(in) :: dPotential1(:,:,:,:)

    !> 1st derivatives of external potential (mOrb, mOrb, nAtom, nSpin)
    real(dp), intent(in) :: dPotential2(:,:,:,:)

    !> Mixed 2nd derivative of external potential (mOrb, mOrb, nAtom, nSpin)
    real(dp), intent(in), allocatable :: d2Potential(:,:,:,:)

    !> Unperturbed eigenvalues of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

    !> Number of (partly) filled states in each [nIndepHam,kpt]
    integer, intent(in) :: nFilled(:,:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Sparse hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Chemical species
    integer, intent(in) :: species(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Map from image atoms to the original unique atom in the central cell
    integer, intent(in) :: img2CentCell(:)

    real(dp), allocatable :: work(:,:), work2(:,:), vSparse(:, :)
    real(dp) :: eTerm
    integer :: iOrb, nOrb, nSpin

    nOrb = orb%nOrb
    nSpin = size(ham, dim=2)

    allocate(work(nOrb, nOrb), source=0.0_dp)
    allocate(work2(nOrb, nOrb), source=0.0_dp)

    call unpackHS(work, ham(:,1), neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call adjointLowerTriangle(work)

    work2(:,:) = 0.0_dp
    do iOrb = 1, nFilled(1,1)
      work2(:, iOrb) = dPsiReal2(:,iOrb) * filling(iOrb, 1, 1)**2
    end do

    ! <Psi1 | H | Psi2>
    eTerm = sum(dPsiReal1(:,:nFilled(1,1)) * matmul(work, work2(:,:nFilled(1,1))))
    !write(*,*)eTerm
    d2E = eTerm

    work(:,:) = 0.0_dp
    call unpackHS(work, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call adjointLowerTriangle(work)

    work2(:,:) = 0.0_dp
    do iOrb = 1, nFilled(1,1)
      work2(:, iOrb) = dPsiReal2(:,iOrb) * eigVals(iOrb, 1, 1) * filling(iOrb, 1, 1)**2
    end do

    ! <Psi1 | - e S| Psi2>
    eTerm = -sum(dPsiReal1(:,:nFilled(1,1)) * matmul(work, work2(:,:nFilled(1,1))))
    !write(*,*)eTerm
    d2E = d2E + eTerm

    deallocate(work2)
    allocate(vsparse(size(over), nSpin))

    vSparse(:,:) = 0.0_dp
    call addShift(env, vSparse, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
        & iSparseStart, nAtom, img2CentCell, dPotential1, isInputZero=.true.)

    work(:,:) = 0.0_dp
    call unpackHS(work, vSparse(:,1), neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call adjointLowerTriangle(work)

    ! <Psi0 | v1 | Psi2>
    eTerm = sum(eigVecsReal(:,:nFilled(1,1), 1) * matmul(work, dPsiReal2(:,:nFilled(1,1))))
    write(*,*)eTerm
    d2E = d2E + eTerm

    vSparse(:,:) = 0.0_dp
    call addShift(env, vSparse, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
        & iSparseStart, nAtom, img2CentCell, dPotential2, isInputZero=.true.)

    work(:,:) = 0.0_dp
    call unpackHS(work, vSparse(:,1), neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call adjointLowerTriangle(work)

    ! <Psi1 | v2 | Psi0>
    eTerm = sum(dPsiReal1(:,:nFilled(1,1)) * matmul(work, eigVecsReal(:,:nFilled(1,1), 1)))
    write(*,*)eTerm
    d2E = d2E + eTerm

    !if (allocated(d2Potential)) then

    !  call addShift(env, vSparse, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
    !      & iSparseStart, nAtom, img2CentCell, d2Potential, isInputZero=.true.)

    !  call unpackHS(work, vSparse(:,1), neighbourList%iNeighbour, nNeighbourSK,&
    !      & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    !  call adjointLowerTriangle(work)

    !  ! <Psi0 | v12| Psi0>
    !  eTerm = sum(eigVecsReal(:,:nFilled(1,1), 1)*matmul(work, eigVecsReal(:,:nFilled(1,1), 1)))
    !  write(*,*)eTerm
    !  d2E = d2E + eTerm

    !end if

    !write(*,*)"Psi"
    !write(*,*)eigVecsReal(:,:nFilled(1,1), 1)
    !write(*,*)"Psi'1"
    !write(*,*)dPsiReal1(:,:nFilled(1,1))
    !write(*,*)"Psi'2"
    !write(*,*)dPsiReal2(:,:nFilled(1,1))
    !write(*,*)"V1"
    !write(*,*)dPotential1
    !write(*,*)"V2"
    !write(*,*)dPotential2

  end subroutine secondEDerivative_real

end module dftbp_derivs_energyderivatives
