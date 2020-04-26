!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains geometrical boundary condition information on the calculation
module dftbp_boundaryconditions
  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
#:if WITH_TRANSPORT
  use libnegf_vars, only : TTransPar
#:endif
  implicit none

  private
  public :: TBoundaryConditions, boundaryTypes, BoundaryConditions_init, zAxis

  type :: TBoundaryTypesEnum

    !> real space cluster/molecular
    integer :: cluster = 0

    !> periodic 3D
    integer :: periodic3D = 1

    !> periodic 1D
    integer :: periodic1D = 2

    !> objective single helical operation (effectively periodic1D + twisting)
    integer :: helical = 3

    !> objective two helical operations
    integer :: helical2Op = 4

    !> Rotational symmetry in xy plane of order N
    integer :: rotational = 5

  end type TBoundaryTypesEnum


  !> Container for enumerated geometric boundary types
  type(TBoundaryTypesEnum), parameter :: boundaryTypes = TBoundaryTypesEnum()


  !> Type for geometric boundary conditions
  type :: TBoundaryConditions

    private

    !> Geometry type of the system
    integer :: boundaryType

    !> Lattice vectors for periodic cases
    real(dp), allocatable :: latVec(:,:)

    !> Open boundary conditions for transport or Dyson-type embedding
    logical :: openBoundary

  end type TBoundaryConditions

  !> z direction vector for rotation
  real(dp), parameter :: zAxis(3) = (/0.0_dp,0.0_dp,1.0_dp/)

contains


  !> Initialise the type of boundary condition on the geometry
#:if WITH_TRANSPORT
  subroutine boundaryConditions_init(this, transport, latVec)
#:else
    subroutine boundaryConditions_init(this, latVec)
#:endif

    !> Self case
    type(TBoundaryConditions), intent(inout) :: this

  #:if WITH_TRANSPORT
    !> Transport input parameters
    type(TTranspar), intent(in) :: transport
  #:endif

    !> Lattice vectors for the system
    real(dp), intent(in), optional :: latVec(:,:)

    if (allocated(this%latVec)) then
      deallocate(this%latVec)
    end if

    if (present(latVec)) then

      if (all(shape(latVec) == [3,3])) then

        this%boundaryType = boundaryTypes%periodic3D

      else if (all(shape(latVec) == [3,1])) then

        this%boundaryType = boundaryTypes%helical2Op

      else if (all(shape(latVec) == [2,1])) then

        this%boundaryType = boundaryTypes%helical

      else

        call error("Unknown boundary conditions")

      end if

      this%latVec = latVec

    else

      this%boundaryType = boundaryTypes%cluster

    end if

  #:if WITH_TRANSPORT
    this%openBoundary = transport%defined
  #:else
    this%openBoundary = .false.
  #:endif

  end subroutine boundaryConditions_init


end module dftbp_boundaryconditions
