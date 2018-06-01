!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exporting the functionality we use from the library dftd3.
module ddpcm
  use accuracy
#:if WITH_DDPCM
  use omp_lib
  use ddcosmo, only : ddinit, memfree, ccav, ncav, nylm, eps, eta, iprint, nproc, lmax,&
      & ngrid, iconv, igrad

  use constants, only : pi, Bohr__AA
  use scc
  use dispuffdata
  use message
  use environment
#:endif
  implicit none

  private

  public :: TddPCM
#:if WITH_DDPCM
  public :: ddPCM_init, buildPCMfields, resetPCM
#:endif

  type :: TddPCM

#:if WITH_DDPCM

    !> dielectric constant
    real(dp) :: dielectric

    real(dp) :: scaleFactor

    real(dp) :: regularisation

    !> van der Waals radii for the species present
    real(dp), allocatable :: vdwRadii(:)

    integer :: nCav

    integer :: nylm

    real(dp), allocatable :: phi(:)

    real(dp), allocatable :: psi(:,:)

    real(dp), allocatable :: sigma(:,:)

  contains

    procedure :: buildPCMfields
    procedure :: resetPCM

#:endif

  end type TddPCM

#:if WITH_DDPCM

  !> Whether code was built with ddPCM support
  logical, parameter, public :: withddPCM = .true.

#:else

  !> Whether code was built with ddPCM support
  logical, parameter, public :: withddPCM = .false.

#:endif

contains

#:if WITH_DDPCM

  subroutine ddPCM_init(this, species0, speciesNames)
    type(TddPCM), intent(inout) :: this

    integer, intent(in) :: species0(:)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    integer :: iSp, iAt, nAtom
    logical :: tFound
    real(dp) :: tmp

    nAtom = size(species0)

    if (.not.allocated(this%vdwRadii)) then
      allocate(this%vdwRadii(nAtom))

      do iAt = 1, nAtom
        iSp = species0(iAt)
        call getUffValues(speciesNames(iSp), this%vdwRadii(iAt), tmp, found=tFound)

        if (.not. tFound) then
          call error("Missing van der Waals radius for this species : " // trim(speciesNames(iSp)) )
        end if

      end do
      this%vdwRadii = this%scaleFactor*this%vdwRadii
      !write(*,*)'Radii :',this%vdwRadii * Bohr__AA
    end if

  end subroutine ddPCM_init

  subroutine resetPCM(this, coord0)
    class(TddPCM), intent(inout) :: this

    real(dp), intent(in) :: coord0(:,:)
        integer :: nAtom

    if (allocated(this%phi)) then
      call memfree()
      deallocate(this%phi)
      deallocate(this%psi)
      deallocate(this%sigma)
    end if

    nAtom = size(coord0,dim=2)

    eps = this%dielectric
    eta = this%regularisation
  #:if DEBUG == 0
    iprint = 0
  #:else
    iprint = 2
  #:endif
    nproc = omp_get_max_threads()
    lmax = 6
    ngrid = 110
    iconv = 7
    igrad = 0
    call ddinit(nAtom, coord0(1,:), coord0(2,:), coord0(3,:), this%vdwRadii(:))

    this%nCav = ncav
    this%nylm = nylm
    allocate(this%phi(ncav))
    allocate(this%psi(nylm,nAtom))
    allocate(this%sigma(nylm,nAtom))
    this%phi = 0.0_dp
    this%psi = 0.0_dp
    this%sigma = 0.0_dp

  end subroutine resetPCM

  subroutine buildPCMfields(this, sccCalc, env, coord0, deltaQAtom, atomPotential, eSolvation)

    class(TddPCM), intent(inout) :: this

    !> Module variables for SCC
    type(TScc), intent(inout) :: sccCalc

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    real(dp), intent(in) :: coord0(:,:)

    real(dp), intent(in) :: deltaQAtom(:)

    real(dp), intent(inout) :: atomPotential(:)

    real(dp), intent(out), optional :: eSolvation

    real(dp) :: fac, esolv, xx(1)
    real(dp), parameter :: tokcal=627.509469_dp

    call sccCalc%getInternalElStatPotential(this%phi, env, ccav)

    ! monopole prefactor
    fac = 2.0_dp * sqrt(pi)

    this%psi = 0.0_dp
    this%psi(1,:) = fac * deltaQAtom(:)
    esolv = 0.0_dp

    call cosmo(.false., .true., this%phi, xx, this%psi, this%sigma, esolv)
    atomPotential(:) = atomPotential(:) + this%sigma(1,:)
    if (present(eSolvation)) then
      eSolvation = esolv
    end if

  end subroutine buildPCMfields

#:endif

end module ddpcm
