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
  use dftbp_geometry_boundarycond, only : boundaryCondsEnum, TBoundaryConds
  use dftbp_math_binarysearch, only : isPresent, first, last
  use dftbp_math_sorting, only : merge_multikey
  implicit none

  private
  public :: TNeighborFinder, TNeighborFinder_init, update


  !> ADT for finding neighbours with O(N ln N) algorithm
  type :: TNeighborFinder

    !> Neighbour cutoff radius
    real(dp) :: cutoff = 0.0_dp

    !> Sorting index for atoms
    integer, allocatable :: sortingIndex(:)

    !> Integer key for which box each atom belongs to
    integer, allocatable :: boxIndex(:,:)

    !> Order and which keys to sort on
    integer :: sortFields(3) = [1,2,3]

  contains

    procedure :: updateIndexing
    procedure :: neighbours

  end type TNeighborFinder

contains


  !> Initialise instance
  subroutine TNeighborFinder_init(this, cutoff, nAtoms)

    type(TNeighborFinder), intent(out) :: this

    real(dp), intent(in) :: cutoff

    integer, intent(in) :: nAtoms

    this%cutoff = cutoff

    allocate(this%sortingIndex(nAtoms), source=0)
    allocate(this%boxIndex(3,nAtoms), source=0)

  end subroutine TNeighborFinder_init


  !> Update indexing for atom structure
  subroutine updateIndexing(this, coords, boundaryConds, errStatus)

    !> Instance
    class(TNeighborFinder), intent(inout) :: this

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Boundary conditions on the calculation
    type(TBoundaryConds), intent(in) :: boundaryConds

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    real(dp) :: invCutoff
    integer :: iAt, nAt
    integer, allocatable :: tmp(:,:)

    nAt = size(coords, dim=2)

    @:ASSERT(nAt == size(this%sortingIndex))
    @:ASSERT(all(shape(this%boxIndex) == [3,nAt]))

    invCutoff = 1.0_dp / this%cutoff
    if (any(abs(coords * invCutoff) > real(huge(1) - 1, dp))) then
      @:RAISE_ERROR(errStatus, -1, "Internal error in dftbp_geometry_neighbours:update insufficent&
          & integer model for keys")
    end if
    this%boxIndex(:, :) = nint(coords * invCutoff)

    call merge_multikey(this%sortingIndex, this%boxIndex, this%sortFields)
    allocate(tmp, mold=this%boxIndex)
    tmp = this%boxIndex(:, this%sortingIndex)
    call move_alloc(tmp, this%boxIndex)

  end subroutine updateIndexing


  !> Find the neighbours of atoms
  subroutine neighbours(this, iFirst, iLast,  coords, boundaryConds, errStatus)

    !> Instance
    class(TNeighborFinder), intent(in) :: this

    !> First atom in range to find
    integer, intent(in) :: iFirst

    !> Last atom in range to find
    integer, intent(in) :: iLast

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Boundary conditions on the calculation
    type(TBoundaryConds), intent(in) :: boundaryConds

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer :: iAt, iBox, ix, iy, iz, jLower, jUpper, key(3), nNeighbouringBoxes
    integer :: neighbouringBox(2,27), oldKey(3),tmpKey(3), nAt
    real(dp) :: cutoff2

    cutoff2 = this%cutoff**2
    nAt = size(this%sortingIndex)
    oldKey = huge(1)
    do iAt = iFirst, iLast
      key(:) = this%boxIndex(:, iAt)
      if (any(key /= oldKey)) then
        ! update which boxes to check for neighbouring atoms
        nNeighbouringBoxes = 0
        do ix = -1, 1
          do iy = -1, 1
            do iz = -1, 1
              tmpKey(:) = key + [ix,iy,iz]
              if (isPresent(1, nAt, this%boxIndex, tmpKey, this%sortFields)) then
                nNeighbouringBoxes = nNeighbouringBoxes + 1
                write(*,*)'Found neighbour box'
                call rangeInBox(jLower, jUpper, tmpKey, this%boxIndex, this%sortFields, errStatus)
                @:PROPAGATE_ERROR(errStatus)
                neighbouringBox(:, nNeighbouringBoxes) = [jLower, jUpper]
              end if
            end do
          end do
        end do
        oldKey(:) = key
      end if

      ! Actual neighbour checks
      do iBox = 1, nNeighbouringBoxes

        jLower = neighbouringBox(1, iBox)
        jUpper = neighbouringBox(2, iBox)
        call findPairs(iAt, iAt, jLower, jUpper, coords, this%sortingIndex, cutoff2)

      end do

    end do

  end subroutine neighbours


  !> Update neighbour data
  subroutine update(coords, cutoff, errStatus)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Cutoff distance to generate neighbours
    real(dp), intent(in) :: cutoff

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer, allocatable :: indx(:), atomBoxes(:,:), sortedAtomBoxes(:,:)
    integer :: iAt, nAt, ix, iy, iz, key(3), tmpKey(3), iLower, iUpper, jLower, jUpper
    real(dp) :: invCutoff, cutoff2
    integer, parameter :: fields(3) = [1,2,3]

    invCutoff = 1.0_dp / cutoff
    cutoff2 = cutoff**2

    nAt = size(coords, dim=2)
    allocate(indx(nAt))
    allocate(atomBoxes(3, nAt))

    if (any(abs(coords * invCutoff) > real(huge(1) - 1, dp))) then
      @:RAISE_ERROR(errStatus, -1, "Internal error in dftbp_geometry_neighbours:update insufficent&
          & integer model for keys")
    end if
    atomBoxes(:, :) = nint(coords * invCutoff)

    call merge_multikey(indx, atomBoxes, fields)
    sortedAtomBoxes = atomBoxes(:,indx)
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
              call findPairs(iLower, iUpper, jLower, jUpper, coords, indx, cutoff2)
            end if
          end do
        end do
      end do

      iAt = iUpper + 1

    end do

  end subroutine update


  !> Calculate range of atoms in the box with the relevant key
  subroutine rangeInBox(iLower, iUpper, key, keys, fields, errStatus)

    !> First in the range
    integer, intent(out) :: iLower

    !> Last in the range
    integer, intent(out) :: iUpper

    !> Index key for the box to check
    integer, intent(in) :: key(:)

    !> List of keys for all the atoms
    integer, intent(in) :: keys(:,:)

    !> Sort order fields
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


  !> Find atoms within the cutoff distance from the supplied range of atoms
  subroutine findPairs(iStart, iEnd, jStart, jEnd, coords, indx, cutoff2)

    !> First of the atoms to find neighbours for (typically within one box)
    integer, intent(in) :: iStart

    !> Last of the atoms to find neighbours for (typically within one box)
    integer, intent(in) :: iEnd

    !> First of the atoms to check for neighbours status (typically within another box)
    integer, intent(in) :: jStart

    !> Last of the atoms to check for neighbours status (typically within the other box)
    integer, intent(in) :: jEnd

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:)

    !> Convert from sorted index to the atom numbers
    integer, intent(in) :: indx(:)

    !> Square of the cutoff distance
    real(dp), intent(in) :: cutoff2

    integer :: ii, jj, iAt, jAt
    real(dp) :: r2

    do ii = iStart, iEnd
      iAt = indx(ii)
      do jj = jStart, jEnd
        jAt = indx(jj)
        r2 = sum((coords(:,iAt) - coords(:,jAt))**2)
        if (r2 <= cutoff2) then
          write(*,"(A,2I4,F12.4)")'Neighbours', iAt, jAt, sqrt(r2) * Bohr__AA
        end if
      end do
    end do

  end subroutine findPairs

end module dftbp_geometry_neighbours
