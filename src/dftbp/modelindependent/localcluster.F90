!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains routines for producing local environments around atoms and convert between neighbour map
!> indexed information and local geometries
module dftbp_modelindependent_localcluster
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_math_simplegeometry, only : dist2Segment
  implicit none

  private
  public :: TLocalCluster, TLocalCluster_init

  type TLocalCluster


  end type TLocalCluster

contains

  !> Initialise local clusters of atoms around
  subroutine TLocalCluster_init(this, cutoffs, coords, boundCond, nNeighbour, neighbourList)

    type(TLocalCluster), intent(out) :: this

    real(dp), intent(in) :: cutoffs(:)

    real(dp), intent(in) :: coords(:,:)

    type(TBoundaryConditions), intent(in) :: boundCond

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, allocatable :: nNeighbour(:)

    !> ADT for neighbour parameters
    type(TNeighbourList), allocatable :: neighbourList


  end subroutine TLocalCluster_init

end module dftbp_modelindependent_localcluster
