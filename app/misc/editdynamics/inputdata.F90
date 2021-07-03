!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains data type representing the input data for editrestart
module editdynamics_inputdata
  use dftbp_common_accuracy, only : dp, lc

  implicit none
  private

  public :: TInputData

  !> container for input data constituents
  type TInputData

    logical :: tInitialized = .false.

    logical, allocatable :: tBinRestartFile(:)

    character(lc), allocatable :: restartFileNames(:)

    !> Masses of atoms in system
    real(dp), allocatable :: atomicMasses(:)

    !> Remove translation motion of centre of mass
    logical, public :: tRemoveTranslate

    !> Remove rotation motion around centre of mass
    logical, public :: tRemoveRotate

  end type TInputData

  type TGeometry
    real(dp), allocatable :: masses(:)
    real(dp), allocatable :: coords(:,:)
    real(dp), allocatable :: velocities(:,:)
    real(dp), allocatable :: oldVelocities(:,:)
  end type TGeometry

end module editdynamics_inputdata
