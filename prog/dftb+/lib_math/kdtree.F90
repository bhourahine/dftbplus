!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> kd-data structure for range searches
module dftbp_kdtree
  use dftbp_accuracy, only : dp
  use dftbp_constants
  use dftbp_sorting
  use dftbp_message
  use dftbp_globalenv, only : stdOut
  implicit  none

  private
  public :: tKDTree, newleaf, insertnodes, rangesearch, printtree, destroytree

  !> Tree node
  type :: tKDTree
    private

    !> Index of coordinate on which the split was performed
    integer :: split

    !> Depth in the tree
    integer :: depth

    !> left descendant node from this one
    type (tKDTree), allocatable :: left

    !> right descendant node from this one
    type (tKDTree), allocatable :: right

  end type tKDTree


  !> Blank string to pad out tree printing for aesthetic reasons
  character (*), parameter :: padding = "                                                          "

contains

  !> Add a new leaf to the tree
  subroutine newleaf(t, depth)

    !> tree instance
    type (tKDTree), intent (out), allocatable :: t

    !> Depth in the tree
    integer, intent (in), optional :: depth

    integer :: thisDepth

    if (present(depth)) then
      thisDepth = depth
    else
      thisDepth = 0
    end if

    allocate(t)

    t % depth = thisDepth

  end subroutine newleaf


  !> Function to insert many points in a reasonably balanced way into a tree. Note there is room to
  !> optimise performance as sorting is on all coordinates at each level. Can instead estimate
  !> median from a sample or use a more complex sorting pattern initially. Also can stop tree with
  !> multiple coordinates in leaf nodes instead of extending to single site at each leaf.
  recursive subroutine insertnodes(t, n, coord, depth)

    !> Tree instance
    type (tKDTree), intent (inout), allocatable :: t

    !> Points in the list of coordinates
    integer, intent (inout) :: n(:)

    !> Coordinates to insert
    real(dp), intent(in)  :: coord(:,:)

    !> Depth in the tree
    integer, intent (in), optional :: depth

    integer :: thisDepth, iDim, nPts, median
    integer, allocatable :: indx(:)

    if (present(depth)) then
      thisDepth = depth
    else
      thisDepth = 0
    end if
    t % depth = thisDepth
    iDim = mod(t % depth, size(coord, dim=1)) + 1
    nPts = size(n)

    select case(nPts)
    case(0)
      call error("Shouldn't be here")
    case(1)

      t % split = n(1)

    case default

      allocate(indx(nPts))
      call index_heap_sort(indx, coord(iDim, n))
      n(:) = n(indx)
      deallocate(indx)
      median = nPts / 2
      if (mod(nPts+1, 2) == 0) then
        median = median + 1
      end if

      ! if the there are other points matching the median coordinate, increase to the end of the
      ! degenerate group, as we are using sorted set instead of comparision against point
      do while(median < nPts)
        if ( abs(coord(iDim, n(median)) - coord(iDim, n(median + 1))) < epsilon(1.0-dp)) then
          median = median + 1
        else
          exit
        end if
      end do

      t % split = n(median)

      if (median > 1) then
        if (.not.allocated(t % left)) then
          allocate(t % left)
        end if
        call insertnodes(t % left, n(:median-1), coord, thisDepth+1)
      end if

      if (median < nPts) then
        if (.not.allocated(t % right)) then
          allocate(t % right)
        end if
        call insertnodes(t % right, n(median+1:nPts), coord, thisDepth+1)
      end if

    end select

  end subroutine insertnodes


  !> Perform an orthogonal range search on the tree around a specified point, returning coordinate
  !> labels within the cutOff
  recursive subroutine rangesearch(t, depth, point, cutOff, coord, nChecked, found, nFound, starter)

    !> Tree instance
    type (tKDTree), intent (in), allocatable :: t

    !> Depth in tree
    integer, intent(in) :: depth

    !> Point around which to compare
    real(dp), intent(in) :: point(:)

    !> Side of box containing point
    real(dp), intent(in) :: cutOff

    !> Coordinates
    real(dp), intent(in)  :: coord(:,:)

    !> number of points being checked
    integer, intent(inout) :: nChecked

    !> List of found points
    integer, intent(inout), allocatable :: found(:)

    !> Number of found points
    integer, intent(inout) :: nFound

    !> starting value to check
    integer, intent(in) :: starter

    integer :: iDim
    integer, allocatable :: tmp(:)

    if (.not. allocated(t)) then
      return
    end if

    nChecked = nChecked + 1

    iDim = mod(t % depth, size(coord, dim=1)) + 1

    if ( sum( (point - coord(:, t % split))**2) <= cutOff**2 ) then
      if (t % split > starter) then
        nFound = nFound + 1
        if (nFound > size(found)) then
          allocate(tmp(2*nFound))
          tmp(:nFound -1) = found(:nFound -1)
          tmp(nFound) = t%split
          call move_alloc(tmp, found)
        end if
      end if
    end if

    if ( point(iDim) - cutOff <= coord(iDim, t % split) ) then

      ! explore left side
      call rangesearch(t%left, depth+1, point, cutOff, coord, nChecked, found, nFound, starter)

    end if

    if ( point(iDim) + cutOff >= coord(iDim, t % split) ) then

      ! explore right side
      call rangesearch(t%right, depth+1, point, cutOff, coord, nChecked, found, nFound, starter)

    end if

  end subroutine rangesearch


  !> Print out the tree if requested
  recursive subroutine printtree(t)

    !> tree node
    type (tKDTree), intent (in), allocatable :: t

    if (allocated(t)) then
      write(stdOut,*) padding(:2*t%depth), 'Node ', t % split, t % depth
      if (allocated(t % left)) then
        call printtree(t % left)
      end if
      if (allocated(t % right)) then
        call printtree(t % right)
      end if
    endif

  end subroutine printtree


  !> Destroy the tree
  recursive subroutine destroytree(t)

    !> Tree node instance
    type (tKDTree), allocatable :: t

    if (.not. allocated(t)) then
      return
    end if

    call destroytree(t % left)
    call destroytree(t % right)
    if (allocated(t % left)) then
      deallocate(t % left)
    end if
    if (allocated(t % right)) then
      deallocate(t % right)
    end if

    deallocate(t)

  end subroutine destroytree

end module dftbp_kdtree
