!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains data type representing the input data for setupgeom
module transporttools_inputdata
  use dftbp_transport_negfvars, only : TTransPar
  use dftbp_type_typegeometry, only : TGeometry
  implicit none
  private

  public :: TInputData
  public :: init, destruct

  !> container for input data constituents
  type TInputData
    logical :: tInitialized = .false.
    type(TGeometry) :: geom
    type(TTransPar) :: transpar
  end type TInputData


  !> Initialise the input data
  interface init
    module procedure InputData_init
  end interface init


  !> destroy input data for variables that do not go out of scope
  interface destruct
    module procedure InputData_destruct
  end interface destruct

contains


  !> Mark data structure as initialised
  subroutine InputData_init(this)
    type(TInputData), intent(out) :: this

    this%tInitialized = .true.

  end subroutine InputData_init


  !> destructor for parts that are not cleaned up when going out of scope
  subroutine InputData_destruct(this)
    type(TInputData), intent(inout) :: this

  end subroutine InputData_destruct

end module transporttools_inputdata
