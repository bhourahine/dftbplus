!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for the atomic neighbour lists
module dftbp_neighbourlists
  use dftbp_assert
  use dftbp_accuracy, only : dp, tolSameDist2
  use dftbp_bisect, only : bisection
  use dftbp_message, only : warning
  implicit none

  private

  public :: TNeighbourList, neighbourList_init
  public :: getNrOfNeighbours, getNrOfNeighboursForAll

  !> Contains essential data for the neighbourlist
  type TNeighbourList

    !> index of neighbour atoms
    integer, allocatable :: iNeighbour(:,:)

    !> nr. of neighbours
    integer, allocatable :: nNeighbour(:)

    !> temporary array for neighbour distances
    real(dp), allocatable :: neighDist2(:,:)

    !> cutoff it was generated for
    real(dp) :: cutoff

    !> initialised data
    logical :: initialized = .false.

  end type TNeighbourList
  
contains

  !> Initializes a neighbourlist instance.
  subroutine neighbourList_init(neighbourList, nAtom, nInitNeighbour)

    !> Neighbourlist data.
    type(TNeighbourList), intent(out) :: neighbourList

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Expected nr. of neighbours per atom.
    integer, intent(in) :: nInitNeighbour

  @:ASSERT(.not. neighbourList%initialized)
  @:ASSERT(nAtom > 0)
  @:ASSERT(nInitNeighbour > 0)

    allocate(neighbourList%nNeighbour(nAtom))
    allocate(neighbourList%iNeighbour(0:nInitNeighbour, nAtom))
    allocate(neighbourList%neighDist2(0:nInitNeighbour, nAtom))

    neighbourList%cutoff = -1.0_dp
    neighbourList%initialized = .true.

  end subroutine neighbourList_init


  !> Returns the nr. of neighbours for a given cutoff for all atoms.
  subroutine getNrOfNeighboursForAll(nNeighbourSK, neigh, cutoff)

    !> Contains the nr. of neighbours for each atom on exit.
    integer, intent(out) :: nNeighbourSK(:)

    !> Initialized neighbourlist
    type(TNeighbourList), intent(in) :: neigh

    !> Maximal neighbour distance to consider.
    real(dp),            intent(in) :: cutoff

    integer :: nAtom, iAtom

    nAtom = size(nNeighbourSK)

    @:ASSERT(size(neigh%iNeighbour, dim=2) == nAtom)
    @:ASSERT(size(neigh%nNeighbour) == nAtom)
    @:ASSERT(maxval(neigh%nNeighbour) <= size(neigh%iNeighbour, dim=1))
    @:ASSERT(all(shape(neigh%neighDist2) == shape(neigh%iNeighbour)))
    @:ASSERT(cutoff >= 0.0_dp)

    ! Get last interacting neighbour for given cutoff
    do iAtom = 1, nAtom
      nNeighbourSK(iAtom) = getNrOfNeighbours(neigh, cutoff, iAtom)
    end do

  end subroutine getNrOfNeighboursForAll


  !> Returns the nr. of neighbours for a given atom.
  function getNrOfNeighbours(neigh, cutoff, iAtom) result(nNeighbour)

    !> Intialised neihgborlist.
    type(TNeighbourList), intent(in) :: neigh

    !> Maximal neighbour distance to consider.
    real(dp),            intent(in) :: cutoff

    !> Index of the atom to get the nr. of neighbours for.
    integer, intent(in) :: iAtom

    !> Nr. of neighbours for the specified atom.
    integer :: nNeighbour

    character(len=100) :: strError
    character(len=*), parameter :: formatErr = "('Cutoff (', E16.6, ') greater then last cutoff ',&
        & '(', E13.6, ') passed to updateNeighbourList!')"

  @:ASSERT(cutoff >= 0.0_dp)
  @:ASSERT(iAtom <= size(neigh%nNeighbour))

    ! Issue warning, if cutoff is bigger than than used for the neighbourlist.
    if (cutoff > neigh%cutoff) then
      write (strError, formatErr) cutoff, neigh%cutoff
      call warning(strError)
    end if

    ! Get last interacting neighbour for given cutoff
    call bisection(nNeighbour, neigh%neighDist2(1:neigh%nNeighbour(iAtom), iAtom), cutoff**2,&
        & tolSameDist2)

  end function getNrOfNeighbours


end module dftbp_neighbourlists
