!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Perturbing fields
module dftbp_derivs_potentialresponse
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init

  private
  public :: TPotentialResponse

! need init inc env, inside scc loop, pass as an array of cases above
! rename as is potential

  !> Generic wrapper of for an external field
  type, abstract :: TPotentialResponse

  contains

    !> Are field cases still left to evaluate?
    procedure(workStillToDo), deferred :: workStillToDo

    !> Advance to next field
    procedure(nextField), deferred :: nextField

    !> Returns field being applied to the system that is charge independent
    procedure(nonSccField), deferred :: nonSccField

    !> Returns field being applied to the system that responds to charges
    procedure(sccField), deferred :: sccField

    !> Post-process field perturbation
    procedure(postProcess), deferred :: postProces

    procedure(completeProcess), deferred :: completeProcess

  end type TPotentialResponse

  abstract interface


    !> Are there perturbing fields still to evaluate?
    function workStillToDo(this, env)
      import : TEnvironment, TPotentialResponse

      !> Instance
      type(TPotentialResponse), intent(in) :: this

      !> Environment settings
      type(TEnvironment), intent(inout) :: env

      !> Are there still perturbations to be performed?
      logical :: workStillToDo

    end function workStillToDo


    !> Move to the next field
    subroutine nextField(this, env)
      import : TEnvironment, TPotentialResponse

      !> Instance
      type(TPotentialResponse), intent(in) :: this

    end subroutine nextField


    subroutine nonSccField(this, env, dPotential)
      import : TEnvironment, TPotentialResponse, TPotentials

      !> Instance
      type(TPotentialResponse), intent(in) :: this

      !> Environment settings
      type(TEnvironment), intent(inout) :: env

      !> Potential contribution
      type(TPotentials), intent(out) :: dPotential

      !> Charges
      real(dp), intent(in) :: dq(:,:,:)

      !> Block charges
      real(dp), intent(in) :: dqBlock(:,:,:,:)

    end subroutine nonSccField


    subroutine sccField(this, env, dPotential, dq, dqBlock)
      import : TEnvironment, TPotentialResponse, TPotentials

      !> Instance
      type(TPotentialResponse), intent(in) :: this

      !> Environment settings
      type(TEnvironment), intent(inout) :: env

      !> Potential contribution
      type(TPotentials), intent(out) :: dPotential

      !> Charges
      real(dp), intent(in) :: dq(:,:,:)

      !> Block charges
      real(dp), intent(in) :: dqBlock(:,:,:,:)

    end subroutine sccField

  end interface

end module dftbp_derivs_potentialresponse
