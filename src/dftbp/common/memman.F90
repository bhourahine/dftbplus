!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains various constants for memory management
module dftbp_common_memman
  use dftbp_common_accuracy, only : dp

  implicit none

  private
  public :: incrmntOfArray, resizeMat


contains


  !> figures out how much larger an array should be to minimize reallocations in future if the array
  !> grows more
  pure function incrmntOfArray(currentSize)
    integer :: incrmntOfArray

    !> current array size
    integer, intent(in) :: currentSize

    incrmntOfArray = currentSize + currentSize  / 2 + 1

  end function incrmntOfArray


  !> Resize matrix
  subroutine resizeMat(mat, newSize)

    !> Matrix to resize
    integer, intent(inout), allocatable :: mat(:)

    !> Lower limit for new size of the matrix
    integer, intent(in) :: newSize

    integer :: oldSize
    integer, allocatable :: work(:)

    oldSize = size(mat)
    if (newSize > oldSize) then
      call move_alloc(mat, work)
      allocate(mat(incrmntOfArray(newSize)), source=0)
      mat(:oldSize) = work
    end if

  end subroutine resizeMat

end module dftbp_common_memman
