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
  use dftbp_dftb_potentials, only : TPotentials
  use dftbp_dftb_shift, only : addShift, totalShift
  use dftbp_dftb_sparse2dense, only : unpackHS
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
  subroutine secondEDerivative_real(env, dPsiReal, dPotential, d2Potential, eigvals, eigVecsReal,&
      & over, ham, nNeighbourSK, neighbourList, species, orb, denseDesc, iSparseStart, nAtom,&
      & img2CentCell, work)

    !> Computational environment settings
    type(TEnvironment), intent(inout) :: env

    !> Derivatives of wavefunctions with respect to perturbations (nOrb, nOrb, nPerturb)
    real(dp), intent(in) :: dPsiReal(:,:,:)

    !> 1st derivatives of external potential (nPerturb)
    type(TPotentials), intent(in) :: dPotential(:)

    !> Mixed 2nd derivative of external potential
    type(TPotentials), intent(in), allocatable :: d2Potential

    !> Unperturbed eigenvalues of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

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

    !> Square work matrix
    real(dp), intent(in) :: work(:,:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Map from image atoms to the original unique atom in the central cell
    integer, intent(in) :: img2CentCell(:)

    real(dp), allocatable :: dHamSqr(:,:), dOverSqr(:,:), vSparse(:, :)
    real(dp), allocatable :: eCiReal(:,:,:)

    integer :: nOrb, nSpin

    @:ASSERT(size(dPotential) == 2)

    nOrb = size(work,dim=1)
    nSpin = size(ham, dim=2)

    allocate(dHamSqr(nOrb, nOrb))
    allocate(dOverSqr(nOrb, nOrb))
    allocate(vSparse(size(over), nSpin))

    call unpackHS(dHamSqr, ham(:,1), neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call unpackHS(dOverSqr, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)



    call addShift(env, vSparse, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
        & iSparseStart, nAtom, img2CentCell, dPotential(1)%extBlock, isInputZero=.true.)

    call addShift(env, vSparse, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
        & iSparseStart, nAtom, img2CentCell, dPotential(2)%extBlock, isInputZero=.true.)

    if (allocated(d2Potential)) then
      call addShift(env, vSparse, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, d2Potential%extBlock, isInputZero=.true.)

    end if

  end subroutine secondEDerivative_real

end module dftbp_derivs_energyderivatives
