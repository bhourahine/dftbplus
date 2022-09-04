!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Connects to an external energy/hamiltonian model or provides a dummy external model if no
!> external library is linked
module dftbp_externalmodel
  use iso_c_binding, only : c_int, c_char, c_bool, c_null_char, c_ptr, c_double, c_loc, c_f_pointer
  use dftbp_io_clang, only : fortranChar
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  implicit none

  private
  public :: TExtModelProvides, TExtModel, getExtModelCapabilities, externalModel_init

  !> Return type for capabilities of external electronic structure/energy model
  type TExtModelProvides

    !> Name of external model
    character(:), allocatable :: modelName

    !> Is a hamiltonian provided
    logical :: gives_ham = .false.

    !> Is an overlap provdided, or should a diagonal matrix be substituted
    logical :: gives_overlap = .false.

    !> Are double-counting energies provided (i.e. non-bandstructure)
    logical :: gives_dc_energy = .false.

    !> Is the model self-consistently evaluated
    logical :: is_self_consistent = .false.

    !> Are results for a subset of the atoms provided (allowing MPI parallism on the DFTB+ side)
    logical :: gives_atomic_subset = .false.

    !> Does the model require a communicator?
    logical :: is_mpi = .false.

    !> Internal model to cover parts not included in the external model (also dictates if total
    !> energies/forces/etc. can be calcualted)
    integer :: intModel = hamiltonianTypes%none

    !> Is an internal model used for repulsive/double counting?
    logical :: is_dc_internal = .false.

  end type TExtModelProvides


  !> Type for instance of model
  type TExtModel

    type(TExtModelProvides) :: capabilities

    real(dp) :: cutoff
    integer, allocatable :: nshells(:)

    !> C pointer to internal state of the model (assume that model manages this)
    type(c_ptr) :: modelState

  end type TExtModel

  !> ctype for receiving external capabilities
  type, bind(C) :: extModelAbilities

    logical(c_bool) :: gives_ham
    logical(c_bool) :: gives_overlap
    logical(c_bool) :: gives_dc_energy
    logical(c_bool) :: requires_self_consistency
    logical(c_bool) :: gives_atom_subset
    logical(c_bool) :: requires_mpi

  end type extModelAbilities


  !> C code interface
  interface

    !> External model declared capabilities
    subroutine model_provides(modelname, capabilities)&
        & bind(C, name='dftbp_provided_with')
      import :: extModelAbilities
      import :: c_char
      implicit none
      character(c_char), intent(out) :: modelname(*)
      type(extModelAbilities), intent(out) :: capabilities
    end subroutine model_provides

    !> Initialise external model for calculation
    function model_init(nspecies, speciesnames, maxCutoff, nShells, shells, modelstate, errString)&
        & bind(C, name='initialise_model_for_dftbp')
      import :: c_int, c_ptr, c_char, c_double
      implicit none
      integer(c_int), intent(in) :: nspecies
      type(c_ptr), target, intent(in) :: speciesnames(nspecies)
      real(c_double), intent(out) :: maxCutoff
      integer(c_int), intent(out) :: nShells(:)
      type(c_ptr), target, intent(out) :: shells(*)
      type(c_ptr), target, intent(out) :: modelstate
      character(c_char), intent(out) :: errString(*)
      integer(c_int) :: model_init
    end function model_init

    !> Initialise external model for calculation
    function model_update(modelstate, errString)&
        & bind(C, name='update_model_for_dftbp')
      import :: c_int, c_ptr, c_char
      implicit none
      !> Internal state of model, opaque to us
      type(c_ptr), target, intent(in) :: modelstate
      !> Any returned error string
      character(c_char), intent(out) :: errString(*)
      !> Model state after operation
      integer(c_int) :: model_update
    end function model_update

    !> Clean up memory attached to a C pointer
    subroutine c_free(ptr) bind(c, name="free")
      import :: c_ptr
      implicit none
      !> Pointer to nullify
      type(c_ptr), value :: ptr
    end subroutine c_free

  end interface

  !> Buffer size on the Fortran side
  integer, parameter :: nBufChar = 256

contains

  !> What are the capabilities of the attached external model
  subroutine getExtModelCapabilities(modelProvides, status)

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Capabilities of externally provided hamiltonian/energy model
    type(TExtModelProvides), intent(out) :: modelProvides

    !> Structure for returned model capabilities
    type(extModelAbilities) :: capabilities

    !> Buffer on Fortran side, expecting a null termination somewhere inside of this, or throws an
    !> error
    character(nBufChar) :: modelname = " "

    call model_provides(modelname, capabilities)
    if (.not.isCStringOK(modelname, nBufChar)) then
      @:RAISE_ERROR(status, -1, "External model name string does not fit in buffer")
    end if
    modelProvides%modelName = trim(modelname)

    modelProvides%gives_ham = capabilities%gives_ham
    modelProvides%gives_overlap = capabilities%gives_overlap
    modelProvides%gives_dc_energy = capabilities%gives_dc_energy
    modelProvides%is_self_consistent = capabilities%requires_self_consistency
    modelProvides%gives_atomic_subset = capabilities%gives_atom_subset
    modelProvides%is_mpi = capabilities%requires_mpi

  end subroutine getExtModelCapabilities


  !> Initialise external model for current calculation
  subroutine externalModel_init(extModel, speciesNames, status)

    !> Instance of external model
    type(TExtModel), intent(out) :: extModel

    !> labels of atomic species
    character(mc), intent(in) :: speciesNames(:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    integer :: iStatus, ii
    character(nBufChar) :: errString = " "
    integer(c_int) :: nspecies
    integer(c_int), allocatable :: nShells(:)
    type(c_ptr), dimension(size(speciesnames)) :: speciesPtrs
    character(mc+1), allocatable, target :: stringArray(:)
    real(c_double) :: maxCutoff
    type(c_ptr) :: cptr_shells
    integer, pointer :: fptr_shells(:)

    nspecies = size(speciesNames)
    allocate(stringArray(nspecies))
    allocate(nShells(nspecies))
    allocate(extModel%nShells(nspecies))

    do ii = 1, nspecies
      stringArray(ii) = trim(speciesNames(ii)) // c_null_char
      speciesPtrs(ii) = c_loc(stringArray(ii))
    end do

    iStatus = model_init(nspecies, speciesPtrs, maxCutoff, nShells, cptr_shells,&
        & extModel%modelState, errString)

    if (iStatus /= 0) then
      if (.not.isCStringOK(errString, nBufChar)) then
        @:RAISE_ERROR(status, -1, "External model name string does not fit in buffer")
      end if
      @:RAISE_ERROR(status, iStatus, trim(errString))
    end if

    extModel%cutoff = maxCutoff
    extModel%nShells = nShells

    call c_f_pointer(cptr_shells, fptr_shells, [sum(nShells)])
    write(*,*)'Ang shells:', fptr_shells

    ! clean up
    fptr_shells => null()
    call c_free(cptr_shells)

    iStatus = model_update(extModel%modelState, errString)
    if (iStatus /= 0) then
      if (.not.isCStringOK(errString, nBufChar)) then
        @:RAISE_ERROR(status, -1, "External model error string does not fit in buffer")
      end if
      @:RAISE_ERROR(status, iStatus, trim(errString))
    end if

  end subroutine externalModel_init


  !> Checks if string has a null termination within the expected length
  function isCStringOK(string, nBufferChar)

    !> String to check
    character(c_char), intent(in) :: string(*)

    !> Expected max length of string
    integer, intent(in) :: nBufferChar

    !> Is the string within the length and null terminated
    logical isCStringOK

    integer :: ii

    isCStringOK = .false.
    lpBufCheck: do ii = 1, nBufferChar
      if (string(ii) == c_null_char) then
        isCStringOK = .true.
        exit lpBufCheck
      end if
    end do lpBufCheck

  end function isCStringOK

end module dftbp_externalmodel
