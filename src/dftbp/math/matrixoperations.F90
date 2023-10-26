!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Simple matrix operations for which LAPACK does not have a direct call
module dftbp_math_matrixoperations
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: triangleCopySquareMatrix, blockTriangleCopyHS

  !> Copy lower triangle into the upper triangle of a square matrix, obeying hermitian symmetry if
  !> appropriate
  interface triangleCopySquareMatrix
    module procedure symmetrizeSquareMatrix
    module procedure hermitianSquareMatrix
  end interface triangleCopySquareMatrix

  !> Copy lower triangle into the upper triangle of a square matrix, obeying hermitian symmetry if
  !> appropriate. Atomic on-site blocks are ignored.
  interface blockTriangleCopyHS
    module procedure blockSymmetrizeHS
    module procedure blockHermitianHS
  end interface blockTriangleCopyHS

contains


  !> Hermitian symmetrize a square matrix, leaving the on-site atomic blocks alone, by copying other
  !! triangle
  subroutine blockHermitianHS(square, iAtomStart)

    !> Square form matrix.
    complex(dp), intent(inout) :: square(:, :)

    !> Returns the offset array for each atom.
    integer, intent(in) :: iAtomStart(:)

    integer :: nAtom, iAtom, iStart, iEnd, mOrb

    nAtom = size(iAtomStart, dim=1) - 1
    mOrb = iAtomStart(nAtom+1) - 1

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT((size(square, dim=1) == mOrb) .or. (size(square, dim=1) == 2*mOrb))
    @:ASSERT(size(square, dim=1) == mOrb .or. size(square, dim=1) == 2*mOrb)

    do iAtom = 1, nAtom
      iStart = iAtomStart(iAtom)
      iEnd = iAtomStart(iAtom+1) - 1
      square(iStart:iEnd, iEnd+1:mOrb) = transpose(conjg(square(iEnd+1:mOrb, iStart:iEnd)))
    end do

    if (size(square, dim=1) == 2*mOrb) then
      ! 2 component matrix, symmetrize L shaped part around the block above
      do iAtom = 1, nAtom
        iStart = iAtomStart(iAtom) + mOrb
        iEnd = iAtomStart(iAtom+1) + mOrb - 1
        square(:iStart-1,iStart:iEnd) = transpose(conjg(square(iStart:iEnd,:iStart-1)))
      end do
    end if

  end subroutine blockHermitianHS


  !> Symmetrize a square matrix, leaving the on-site atomic blocks alone, by copying the other
  !! triangle
  subroutine blockSymmetrizeHS(square, iAtomStart)

    !> Square form matrix.
    real(dp), intent(inout) :: square(:, :)

    !> Returns the offset array for each atom.
    integer, intent(in) :: iAtomStart(:)

    integer :: nAtom, iAtom, iStart, iEnd

    nAtom = size(iAtomStart) - 1

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))

    do iAtom = 1, nAtom
      iStart = iAtomStart(iAtom)
      iEnd = iAtomStart(iAtom+1) - 1
      square(iStart:iEnd, iEnd+1:) = transpose(square(iEnd+1:, iStart:iEnd))
    end do

  end subroutine blockSymmetrizeHS


  !> Copy lower triangle to upper for a square matrix.
  subroutine symmetrizeSquareMatrix(matrix)

    !> matrix to symmetrize
    real(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
    do ii = 1, matSize -1
      matrix(ii, ii+1:) = matrix(ii+1:, ii)
    end do

  end subroutine symmetrizeSquareMatrix


  !> Copy lower triangle to upper for a square matrix with Hermitian symmetry
  subroutine hermitianSquareMatrix(matrix)

    !> matrix to symmetrize
    complex(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
    do ii = 1, matSize -1
      matrix(ii, ii+1:) = conjg(matrix(ii+1:, ii))
    end do

  end subroutine hermitianSquareMatrix

end module dftbp_math_matrixoperations
