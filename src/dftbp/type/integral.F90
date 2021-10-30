!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Data types to handle overlap related integrals
module dftbp_type_integral
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TIntegral, TIntegral_init


  !> Container to store overlap related integrals
  type :: TIntegral

    !> Overlap integrals in atomic block sparse form
    real(dp), allocatable :: overlap(:)

    !> Dipole integrals with operator on ket function
    real(dp), allocatable :: dipoleKet(:, :)

    !> Dipole integrals with operator on bra function
    real(dp), allocatable :: dipoleBra(:, :)

    !> Quadrupole integrals with operator on ket function
    real(dp), allocatable :: quadrupoleKet(:, :)

    !> Quadrupole integrals with operator on bra function
    real(dp), allocatable :: quadrupoleBra(:, :)

    !> Real Hamiltonian integrals in atomic block sparse form
    real(dp), allocatable :: hamiltonian(:, :)

    !> Imaginary Hamiltonian integrals in atomic block sparse form
    real(dp), allocatable :: iHamiltonian(:, :)

    !> Real Hamiltonian integrals in atomic block sparse form transformed by a vector potential
    real(dp), allocatable :: hamiltonianGauged(:,:)

    !> Imaginary Hamiltonian integrals in atomic block sparse form transformed by a vector potential
    real(dp), allocatable :: iHamiltonianGauged(:,:)

    !> Real Overlap integrals in atomic block sparse form transformed by a vector potential (London)
    real(dp), allocatable :: overlapGauged(:)

    !> Imaginary Overlap integrals in atomic block sparse form transformed by a vector potential
    !> (London)
    real(dp), allocatable :: iOverlapGauged(:)

  contains

    !> Reallocates size for the integrals
    procedure :: reallocate

  end type TIntegral


contains


  !> Initializier for integral container
  subroutine TIntegral_init(this, nSpin, isReHam, isImHam, nDipole, nQuadrupole, isGauged)

    !> Instance of the integral container
    type(TIntegral), intent(out) :: this

    !> Number of spins channels in the system
    integer, intent(in) :: nSpin

    !> Allocate space for real Hamiltonian
    logical, intent(in) :: isReHam

    !> Allocate space for imaginary Hamiltonian
    logical, intent(in) :: isImHam

    !> Number of dipole moment components included
    integer, intent(in) :: nDipole

    !> Number of quadrupole moment components included
    integer, intent(in) :: nQuadrupole

    !> Is a gauge field present, leading to complext phases of the basis
    logical, intent(in) :: isGauged

    if (isReHam) then
      allocate(this%hamiltonian(0, nSpin))
      if (isGauged) then
        allocate(this%hamiltonianGauged(0, nSpin))
        allocate(this%iHamiltonianGauged(0, nSpin))
      end if
    end if
    if (isImHam) then
      allocate(this%iHamiltonian(0, nSpin))
    end if
    allocate(this%overlap(0))
    if (isGauged) then
      allocate(this%overlapGauged(0))
      allocate(this%iOverlapGauged(0))
    end if
    if (nDipole > 0) then
      allocate(this%dipoleKet(nDipole, 0))
      allocate(this%dipoleBra(nDipole, 0))
    end if
    if (nQuadrupole > 0) then
      allocate(this%quadrupoleKet(nQuadrupole, 0))
      allocate(this%quadrupoleBra(nQuadrupole, 0))
    end if

  end subroutine TIntegral_init


  !> Re-size storage for integrals
  subroutine reallocate(this, sparseSize, isREKS)

    !> Instance
    class(TIntegral), intent(inout) :: this

    !> Size of the sparse overlap
    integer, intent(in) :: sparseSize

    !> Is this a REKS calculation
    logical, intent(in) :: isREKS

    integer :: nSpin, nDipole, nQuadrupole

    nSpin = size(this%hamiltonian, dim=2)

    if (.not. isREKS) then
      deallocate(this%hamiltonian)
      allocate(this%hamiltonian(sparseSize, nSpin))
    end if

    deallocate(this%overlap)
    allocate(this%overlap(sparseSize))

    if (allocated(this%iHamiltonian)) then
      deallocate(this%iHamiltonian)
      allocate(this%iHamiltonian(sparseSize, nSpin))
    end if

    if (allocated(this%hamiltonianGauged)) then
      deallocate(this%hamiltonianGauged)
      allocate(this%hamiltonianGauged(sparseSize, nSpin))
      deallocate(this%iHamiltonianGauged)
      allocate(this%iHamiltonianGauged(sparseSize, nSpin))
    end if

    if (allocated(this%overlapGauged)) then
      deallocate(this%overlapGauged)
      allocate(this%overlapGauged(sparseSize))
      deallocate(this%iOverlapGauged)
      allocate(this%iOverlapGauged(sparseSize))
    end if

    if (allocated(this%dipoleKet)) then
      nDipole = size(this%dipoleKet, 1)
      deallocate(this%dipoleBra, this%dipoleKet)
      allocate(this%dipoleKet(nDipole, sparseSize))
      allocate(this%dipoleBra(nDipole, sparseSize))
    end if
    if (allocated(this%quadrupoleKet)) then
      nQuadrupole = size(this%quadrupoleKet, 1)
      deallocate(this%quadrupoleBra, this%quadrupoleKet)
      allocate(this%quadrupoleKet(nQuadrupole, sparseSize))
      allocate(this%quadrupoleBra(nQuadrupole, sparseSize))
    end if

  end subroutine reallocate

end module dftbp_type_integral
