!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for molecular boundary conditions
module dftbp_cluster
  use dftbp_accuracy
  use dftbp_sorting
  use dftbp_linkedlist
  use dftbp_constants
  use dftbp_memman
  use dftbp_neighbourlists
  use dftbp_matrixindexing
  implicit none

  private

  public :: clusterNeighbourList

contains


  !> Updates the neighbour list for a given geometry.  The neighbourlist for the given cutoff is
  !> calculated. Arrays are resized if necessary. The neighbour list determination is a simple O(N
  !> log N) algorithm, calculating the distance between the possible atom pairs. Note, can be O(N)
  !> with a better scaling sort algorithm, and eventually the index function will exceed the integer
  !> type for very large ranges of coordinates.
  subroutine clusterNeighbourList(coord, img2CentCell, iCellVec, neigh, nAllAtom,&
      & coord0, cutoff, symmetric, species0,  species)

    !> Coordinates of all interacting atoms on exit
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Mapping on atoms in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Shift vector index for every interacting atom
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> Updated neighbour list.
    type(TNeighbourList), intent(inout) :: neigh

    !> Number of all interacting atoms
    integer, intent(out) :: nAllAtom

    !> Coordinates of the objects in the central cell.
    real(dp), intent(in) :: coord0(:,:)

    !> Cutoff radius for the interactions.
    real(dp), intent(in) :: cutoff

    !> Whether the neighbour list should be symmetric or not (default)
    logical, intent(in), optional :: symmetric

    !> Species of the atoms in the central cell
    integer, intent(in), optional :: species0(:)

    !> Species of all interacting atoms on exit.
    integer, allocatable, intent(inout), optional :: species(:)

    integer :: ii, jj, kk, iAt1, iAt2, nAtom, nAtomNear, nn1, maxNeighbour
    type(TListInt), allocatable :: candidates(:)
    integer, allocatable :: iSurroundAt(:), boxIndex(:), indx(:)
    real(dp) ::  minimum(3), maximum(3), diff(3), cutPlus, cTmp(3), rTmp
    integer :: iBox
    logical :: symm

    nAtom = size(coord0, dim=2)

    coord = coord0
    if (present(species)) then
    @:ASSERT(present(species0))
      species = species0
    end if
    do iAt1 = 1, nAtom
      img2CentCell(iAt1) = iAt1
    end do
    iCellVec(:) = 0
    nAllAtom = nAtom

    symm = .false.
    if (present(symmetric)) then
      symm = symmetric
    end if

    neigh%cutoff = cutoff

    ! Clean out arrays
    !  (Every atom is the 0th neighbour of itself with zero distance square.)
    neigh%nNeighbour(:) = 0
    neigh%iNeighbour(:,:) = 0
    do ii = 1, nAtom
      neigh%iNeighbour(0, ii) = ii
    end do
    neigh%neighDist2(:,:) = 0.0_dp

    cutPlus = cutOff * 1.1_dp

    allocate(candidates(nAtom))
    do ii = 1, nAtom
      call init(candidates(ii))
    end do

    minimum(:) = minval(coord0, dim=2)
    maximum(:) = maxval(coord0, dim=2)
    diff(:) = maximum - minimum

    allocate(boxIndex(27*nAtom))
    ! add displaced copies of atom into surrounding boxes
    iAt1 = 1
    do ii = -1, 1
      do jj = -1, 1
        do kk = -1, 1
          do iAt2 = 1, nAtom
            cTmp = coord0(:,iAt2)
            cTmp(1) = cTmp(1) + ii * cutPlus*.5
            cTmp(2) = cTmp(2) + jj * cutPlus*.5
            cTmp(3) = cTmp(3) + kk * cutPlus*.5
            boxIndex(iAt1) = boxIndx(cTmp, minimum, diff, cutPlus)
            iAt1 = iAt1 + 1
          end do
        end do
      end do
    end do
    allocate(indx(27*nAtom))
    call index_heap_sort(indx, boxIndex)

    do ii = 1, 27*nAtom
      iBox = boxIndex(indx(ii))
      iAt1 = mod(indx(ii)-1, nAtom)+1
      do jj = ii-1, 1, -1
        if (iBox == boxIndex(indx(jj))) then
          iAt2 = mod(indx(jj)-1,nAtom)+1
          if (iAt1 /= iAt2) then
            call append(candidates(iAt1), iAt2)
          end if
        else
          exit
        end if
      end do
      do jj = ii+1, nAtom
        if (iBox == boxIndex(indx(jj))) then
          iAt2 = mod(indx(jj)-1,nAtom)+1
          if (iAt1 /= iAt2) then
            call append(candidates(iAt1), iAt2)
          end if
        else
          exit
        end if
      end do
    end do

    deallocate(indx)
    deallocate(boxIndex)

    maxNeighbour = 0
    do iAt1 = 1, nAtom
      if (len(candidates(iAt1)) == 0) then
        cycle
      end if
      allocate(iSurroundAt(len(candidates(iAt1))))
      call asArray(candidates(iAt1),iSurroundAt)
      call heap_sort(iSurroundAt)
      nAtomNear = unique(iSurroundAt)
      maxNeighbour = max(maxNeighbour, nAtomNear)
      if (nAtomNear > size(neigh%iNeighbour, dim=1)) then
        call reallocateArrays3(neigh%iNeighbour, neigh%neighDist2, incrmntOfArray(nAtomNear))
      end if
      do ii = 1, nAtomNear
        iAt2 = iSurroundAt(ii)
        if (.not.symm .and. iAt2 < iAt1) then
          cycle
        end if
        rTmp = sum((coord0(:,iAt1) - coord0(:,iAt2))**2)
        if (rTmp <= cutOff**2) then
          neigh%nNeighbour(iAt1) = neigh%nNeighbour(iAt1) + 1
          neigh%iNeighbour(neigh%nNeighbour(iAt1),iAt1) = iAt2
          neigh%neighDist2(neigh%nNeighbour(iAt1),iAt1) = rTmp
        end if
      end do
      call destruct(candidates(iAt1))
      deallocate(iSurroundAt)
    end do

    ! Sort neighbours of each atom by distance
    allocate(indx(maxNeighbour))
    do iAt1 = 1, nAtom
      nn1 = neigh%nNeighbour(iAt1)
      call index_heap_sort(indx(1:nn1), neigh%neighDist2(1:nn1, iAt1), tolSameDist2)
      neigh%iNeighbour(1:nn1, iAt1) = neigh%iNeighbour(indx(:nn1), iAt1)
      neigh%neighDist2(1:nn1, iAt1) = neigh%neighDist2(indx(:nn1), iAt1)
    end do

  end subroutine clusterNeighbourList

  !> Function to map coordinate into a box number
  pure function boxIndx(coord, minimum, diff, cutOff)

    !> Coordinate
    real(dp), intent(in) :: coord(3)

    !> Minimum of all the coordinates in the system
    real(dp), intent(in) :: minimum(3)

    !> range of the coordinates in the system
    real(dp), intent(in) :: diff(3)

    !> Box side
    real(dp), intent(in) :: cutOff

    !> Return
    integer :: boxIndx

    boxIndx = nint((coord(1)-minimum(1))/cutOff)&
        & + nint((coord(2)-minimum(2))/cutOff)*ceiling((diff(1)+cutoff)/cutOff)&
        & + nint((coord(3)-minimum(3))/cutOff)*ceiling((diff(1)+cutoff)/cutOff)&
        & * ceiling((diff(2)+cutoff)/cutOff)

  end function boxIndx

end module dftbp_cluster
