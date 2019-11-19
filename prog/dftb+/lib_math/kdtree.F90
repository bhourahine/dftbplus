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
  use dftbp_sorting
  use dftbp_message
  use dftbp_globalenv, only : stdOut
  implicit  none

  private
  public :: tKDTree, newleaf, insertnodes, printtree, destroytree

  !> Tree node
  type :: tKDTree
    private
    integer :: split
    integer :: depth
    type (tKDTree), allocatable :: left, right
  end type tKDTree

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


  !> Function to insert many points in a balanced way into the tree
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


  recursive subroutine printtree(t)
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

  recursive subroutine destroytree(t)

    type (tKDTree), allocatable :: t

    if (.not. allocated(t)) then
      return
    end if

    !write(*,*)padding(:2*t%depth), 'Node with ', t % split, t % depth

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


!program main
!  use class_tree
!  implicit none
!
!  integer :: ii, n, nPts
!  type (tKDTree), allocatable :: t
!
!  allocate(t)
!  read(*,*)nPts
!  do ii = 1, nPts
!    read(*,*)n
!    if (ii == 1) then
!      call newleaf(t, n)
!    else
!      call insertnode(t,n)
!    end if
!  end do
!  write(*,*)'Tree'
!  call printtree(t)
!  write(*,*)'Destroy'
!  call destroytree(t)
!
!end program main
