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
  implicit none

  private
  public :: TBoundaryConditions, boundaryTypes, BoundaryConditions_init

  type :: TBoundaryTypesEnum

    !> real space cluster/molecular
    integer :: cluster = 0

    !> periodic 3D
    integer :: periodic3D = 1

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

  contains

    procedure :: periodic
    procedure :: cluster

  end type TBoundaryConditions

contains

  !> Initialise the type of boundary condition on the geometry
  subroutine BoundaryConditions_init(this, latticeVecs)

    type(TBoundaryConditions), intent(out) :: this

    real(dp), intent(in), allocatable :: latticeVecs(:,:)

    if (allocated(latticeVecs)) then
      this%boundaryType = boundaryTypes%periodic3D
      this%latVec = latticeVecs
    else
      this%boundaryType = boundaryTypes%cluster
    end if

  end subroutine BoundaryConditions_init


  !> Is this a 3D periodic geometry
  function periodic(bc)

    !> Instance
    class(TBoundaryConditions), intent(in) :: bc

    !> Return
    logical :: periodic

    periodic = (bc%boundaryType == boundaryTypes%periodic3D)

  end function periodic


  !> Is this a cluster geometry
  function cluster(bc)

    !> Instance
    class(TBoundaryConditions), intent(in) :: bc

    !> Return
    logical :: cluster

    cluster = (bc%boundaryType == boundaryTypes%cluster)

  end function cluster

  
end module dftbp_boundaryconditions
