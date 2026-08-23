!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains routines for neighbour finding
module dftbp_geometry_neighbours
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : Bohr__AA
  use dftbp_common_status, only : TStatus
  use dftbp_math_binarysearch, only : isPresent, first, last
  use dftbp_math_sorting, only : merge_multikey
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


  subroutine update(coords, cutoff, errStatus)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Cutoff distance to generate neighbours
    real(dp), intent(in) :: cutoff

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer, allocatable :: indx(:), atomBoxes(:,:), sortedAtomBoxes(:,:)
    integer :: iAt, nAt, ix, iy, iz, key(3), tmpKey(3), iLower, iUpper, jLower, jUpper
    real(dp) :: invCutoff
    integer, parameter :: fields(3) = [1,2,3]

    invCutoff = 1.0_dp / cutoff

    nAt = size(coords, dim=2)
    allocate(indx(nAt))
    !allocate(invIndx(nAt))
    allocate(atomBoxes(3, nAt))

    if (any(abs(coords * invCutoff) > real(huge(1) - 1, dp))) then
      @:RAISE_ERROR(errStatus, -1, "Internal error in dftbp_geometry_neighbours:update insufficent&
          & integer model for keys")
    end if
    atomBoxes(:, :) = nint(coords * invCutoff)

    call merge_multikey(indx, atomBoxes, fields)
    sortedAtomBoxes = atomBoxes(:,indx)
    !invIndx(indx) = [(iAt, iAt = 1, nAt)]
    write(*,*)'Range in Box'

    iAt = 1
    do while (iAt <= size(indx))

      write(*,*)'Range call', iAt, ':', sortedAtomBoxes(:,iAt)
      if (isPresent(1, size(indx), sortedAtomBoxes, sortedAtomBoxes(:,iAt), fields)) then
        call rangeInBox(iLower, iUpper, sortedAtomBoxes(:,iAt), sortedAtomBoxes, fields, errStatus)
      else
        @:RAISE_ERROR(errStatus, -1, "Internal error in dftbp_geometry_neighbours")
      end if
      write(*,*)'Box atom range ', iLower, iUpper

      ! Find Neighbouring boxes
      key(:) = sortedAtomBoxes(:, iAt)
      write(*,*)'Atom', iAt, ' is in :', key
      do ix = -1, 1
        do iy = -1, 1
          do iz = -1, 1
            !if (all([ix,iy,iz] == 0)) cycle ! could test central box, making logic simpler
            tmpKey(:) = sortedAtomBoxes(:, iAt)
            tmpKey(:) = tmpKey + [ix,iy,iz]
            if (isPresent(1, size(indx), sortedAtomBoxes, tmpKey, fields)) then
              write(*,*)'Found neighbour box'
              call rangeInBox(jLower, jUpper, tmpKey, sortedAtomBoxes, fields, errStatus)
              @:PROPAGATE_ERROR(errStatus)
              write(*,*)'Neighbour atom range ', jLower, jUpper
              call neighbours(iLower, iUpper, jLower, jUpper, coords, indx, cutoff)
            end if
          end do
        end do
      end do

      iAt = iUpper + 1

    end do

  end subroutine update


  subroutine rangeInBox(iLower, iUpper, key, keys, fields, errStatus)

    integer, intent(out) :: iLower

    integer, intent(out) :: iUpper

    !> Index key for the box to check
    integer, intent(in) :: key(:)

    integer, intent(in) :: keys(:,:)

    integer, intent(in) :: fields(:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer :: nn

    nn = size(keys, dim=2)
    iLower = first(1, nn, keys, key, fields)
    if (iLower == -1) then
      @:RAISE_ERROR(errStatus, -1, "Internal error in dftbp_geometry_neighbours:rangeInBox lower&
          & range")
    end if
    iUpper = last(1, nn, keys, key, fields)
    if (iUpper == -1) then
      @:RAISE_ERROR(errStatus, -1, "Internal error in dftbp_geometry_neighbours:rangeInBox upper&
          & range")
    end if

  end subroutine rangeInBox


  subroutine neighbours(iStart, iEnd, jStart, jEnd, coords, indx, cutoff)

    integer, intent(in) :: iStart, iEnd, jStart, jEnd
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: indx(:)
    real(dp), intent(in) :: cutoff

    integer :: ii, jj, iAt, jAt
    real(dp) :: r2, d2

    d2 = cutoff ** 2
    do ii = iStart, iEnd
      iAt = indx(ii)
      do jj = jStart, jEnd
        jAt = indx(jj)
        r2 = sum((coords(:,iAt) - coords(:,jAt))**2)
        if (r2 <= d2) then
          write(*,"(A,2I4,F12.4)")'NeighboursX', iAt, jAt, sqrt(r2) * Bohr__AA
        end if
      end do
    end do

  end subroutine neighbours

end module dftbp_geometry_neighbours
