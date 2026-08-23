!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for neighbour finding
module dftbp_geometry_neighbours
  use dftbp_common_accuracy, only : dp
  use dftbp_math_binarysearch, only : search_int_multikey
  use dftbp_math_sorting, only : merge_multikey, multicompare_int
  implicit none

  private
  public ::TNeighborFinder, update

  type :: TNeighborFinder

    !> Neighbor cutoff radius
    real(dp) :: cutoff = 0.0_dp

    !> Grid dimensions along each fractional axis to check
    integer :: num_cells(3) = 0

  end type TNeighborFinder

contains

  subroutine TNeighborFinder_init()

  end subroutine TNeighborFinder_init


  subroutine update(coords, cutoff)

    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: cutoff

    integer, allocatable :: indx(:), atomBox(:,:), sortedBoxes(:,:)
    integer :: iAt, nAt, iBox, jBox, ix, iy, iz, key(3), tmpKey(3), iLower, iUpper
    real(dp) :: invCutoff

    invCutoff = 1.0_dp / cutoff

    nAt = size(coords, dim=2)
    allocate(indx(nAt))
    allocate(atomBox(3, nAt))

    write(*,*)'Initial boxes (keys)'
    do iAt = 1, nAt
      atomBox(:, iAt) = nint(coords(:,iAt) * invCutoff)
      write(*,*)iAt, atomBox(:, iAt)
    end do

    call merge_multikey(indx, atomBox)
    write(*,*)'Sort keys'
    sortedBoxes = atomBox(:,indx)
    do iAt = 1, nAt
      write(*,*)iAt, indx(iAt), ':', sortedBoxes(:, iAt)
    end do

    write(*,*)'Assignments'
    do iAt = 1, nAt

    end do

    write(*,*)'Range in Box'
    iAt = 1
    do while (iAt <= size(indx))

      !if (multicompare_int(sortedBoxes(:, jBox), key)) exit
      !write(*,*)'Upper', jBox - 1

      write(*,*)'Range call'
      call rangeInBox(iLower, iUpper, iAt, sortedBoxes)
      write(*,*)'Atom range ', iLower, iUpper

      ! Find Neighbouring boxes
      key(:) = sortedBoxes(:, iAt)
      write(*,*)'Key', key
      write(*,*)iAt, ' is in :', sortedBoxes(:, iAt)
      do ix = -1, 1
        do iy = -1, 1
          do iz = -1, 1
            if (all([ix,iy,iz] == 0)) cycle
            tmpKey(:) = sortedBoxes(:, iAt)
            tmpKey(:) = tmpKey + [ix,iy,iz]
            call search_int_multikey(jBox, sortedBoxes, multicompare_int, tmpKey)
            jBox = jBox + 1
            if (testRangeValidity(jBox, size(indx))) then
              if (any(abs(sortedBoxes(:, jBox) - tmpKey) /= 0)) cycle
              write(*,*)'Neighbour box ', jBox, ':', sortedBoxes(:, jBox)
            end if
          end do
        end do
      end do

      iAt = iUpper + 1

      write(*,*)

    end do

  end subroutine update


  subroutine rangeInBox(iLower, iUpper, initBox, sortedKeys)

    integer, intent(out) :: iLower
    integer, intent(out) :: iUpper
    integer, intent(in) :: initBox
    integer, intent(in) :: sortedKeys(:,:)

    integer :: key(3), indxBox, jBox, n

    n = size(sortedKeys, dim=2)

    key(:) = sortedKeys(:, initBox)
    ! redundant to search this if it is already the lower values of the range:
    call search_int_multikey(indxBox, sortedKeys, multicompare_int, key)
    iLower = indxBox + 1
    if (.not. testRangeValidity(iLower, n)) then
      write(*,*)'Internal err'
    end if
    key(3) = key(3) + 1 ! as least significant change is on z index
    call search_int_multikey(indxBox, sortedKeys, multicompare_int, key)
    iUpper = indxBox
    if (.not. testRangeValidity(iUpper, n)) then
      write(*,*)'Internal err'
    end if

  end subroutine rangeInBox


  pure function testRangeValidity(ii, n)

    integer, intent(in) :: ii
    integer, intent(in) :: n

    logical testRangeValidity

    testRangeValidity = ii >= 1 .and. ii <= n

  end function testRangeValidity

end module dftbp_geometry_neighbours
