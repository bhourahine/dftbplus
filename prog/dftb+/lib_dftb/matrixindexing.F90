!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:set ATTRIB = [('1', ':', 'nElem'), ('2', ':, :', 'nElem, nSpin'), ('single', ':', 'nElem')]

!> Contains subroutines for managing neighbour data for large sparse matrices (H, S)
module dftbp_matrixindexing
  use dftbp_accuracy
  use dftbp_commontypes
  implicit none

  private

  public :: reallocateHS, buildSquaredAtomIndex, getSparseDescriptor

  !> resize sparse arrays
  interface reallocateHS
  #:for NAME, _, _ in ATTRIB
    module procedure reallocateHS_${NAME}$
  #:endfor
  end interface reallocateHS

contains

#:for NAME, SIZE, DIMS in ATTRIB
  !> Allocate (reallocate) space for sparse matrices
  subroutine reallocateHS_${NAME}$(ham,&
    #:if NAME != 'single'
      & over,&
    #:endif
      & iPair, iNeighbour, nNeighbourSK, orb, img2Centcell)

    !> Hamiltonian
    real(dp), allocatable, intent(inout):: ham(${SIZE}$)

  #:if NAME != 'single'
    !> Overlap matrix
    real(dp), allocatable, intent(inout) :: over(:)
  #:endif

    !> Pair indexing array (specifying the offset for the interaction between atoms in the central
    !> cell and their neighbours)
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> List of neighbours for each atom in the central cell. (Note: first index runs from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom in the central cell.
    integer, intent(in) :: nNeighbourSK(:)

    !> Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> array mapping images of atoms to originals in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> nr. atoms in the central cell
    integer :: nAtom

    !> nr. of elements in the sparse H/S before and after resizing
    integer :: nOldElem, nElem

    !> nr. of max. possible neighbours (incl. itself)
    integer :: mNeighbour

    integer :: ind, iAt1, iNeigh1, nOrb1

  #:if NAME == '2'
    integer :: nSpin

    nSpin = size(ham, dim=2)
  #:endif

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
  #:if NAME != 'single'
    @:ASSERT(allocated(over))
    @:ASSERT(size(over) == nOldElem)
  #:endif
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbour-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind + nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1, iAt1)))
      end do
    end do
    nElem = ind
    if (nElem > nOldElem) then
      deallocate(ham)
      allocate(ham(${DIMS}$))
      ham(${SIZE}$) = 0.0_dp
    #:if NAME != 'single'
      deallocate(over)
      allocate(over(nElem))
      over(:) = 0.0_dp
    #:endif
    end if

  end subroutine reallocateHS_${NAME}$

#:endfor


  !> Builds an atom offset array for the squared hamiltonain/overlap.
  subroutine buildSquaredAtomIndex(iAtomStart, orb)

    !> Returns the offset array for each atom.
    integer, intent(out) :: iAtomStart(:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer :: ind, iAt1
    integer :: nAtom

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(iAtomStart) == (/ nAtom + 1 /)))

    ind = 1
    do iAt1 = 1, nAtom
      iAtomStart(iAt1) = ind
      ind = ind + orb%nOrbAtom(iAt1)
    end do
    iAtomStart(nAtom+1) = ind

  end subroutine buildSquaredAtomIndex


  !> Calculate indexing array and number of elements in sparse arrays like the real space overlap
  subroutine getSparseDescriptor(iNeighbour, nNeighbourSK, img2CentCell, orb, iPair, sparseSize)

    !> Neighbours of each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing for mapping image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse array indexing for the start of atomic blocks in data structures
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> Total number of elements in a sparse structure (ignoring extra indices like spin)
    integer, intent(out) :: sparseSize

    integer :: nAtom, mNeighbour
    integer :: ind, iAt1, nOrb1, iNeigh1, nOrb2

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)

  @:ASSERT(allocated(iPair))
  @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0 : mNeighbour - 1, nAtom))
      iPair(:,:) = 0
    end if
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        nOrb2 = orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1, iAt1)))
        ind = ind + nOrb1 * nOrb2
      end do
    end do
    sparseSize = ind

  end subroutine getSparseDescriptor

end module dftbp_matrixindexing
