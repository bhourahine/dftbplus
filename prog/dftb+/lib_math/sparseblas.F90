!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> DFTB+ sparse structure operations for a few blas routines
module dftbp_sparseblas
  use dftbp_accuracy, only : dp
  use dftbp_assert
  use dftbp_periodic, only : TNeighbourList
  implicit none

  private
  public :: dftbSYMV, dftbSYMM
  
  !> Routine for symmetric sparse matrix vector multiply
  interface dftbSYMV
     module procedure SYMV_gamma
     ! module procedure HEMV_kpt
  end interface

  !> Multiplication between a symmetric sparse matrix and a general dense matrix
  interface dftbSYMM
    module procedure SYMM_gamma
  end interface dftbSYMM
  
contains

  !> Symmetric sparse matrix multiplied with a vector y = alpha * A * x + beta * y
  subroutine SYMV_gamma(y, A, x, iNeighbour, nNeighbour, img2CentCell, iPair, iSquare, mOrb, alpha,&
      & beta)

    !> Vector
    real(dp), intent(inout) :: y(:)

    !> Symmetric matrix in sparse format
    real(dp), intent(in) :: A(:)

    !> Vector to multiply
    real(dp), intent(in) :: x(:)

    !> Neighbouring atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbouring atoms
    integer, intent(in) :: nNeighbour(:)

    !> Atom indices into central cell atom numbers
    integer, intent(in) :: img2CentCell(:)

    !> Indexing for sparse matrix
    integer, intent(in) :: iPair(0:,:)

    !> Indexing for dense matrinx
    integer, intent(in) :: iSquare(:)

    !> Maximum number of orbitals in an atomic block
    integer, intent(in) :: mOrb

    !> Scaling factor for A x, 1 if not set
    real(dp), optional, intent(in) :: alpha

    !> Scaling factor for incomming y, 0 if not set
    real(dp), optional, intent(in) :: beta
    
    integer :: iOrig, ix, iy, jx, jy, iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, ii
    real(dp) :: sqrTmp(mOrb,mOrb), alphaTmp
    
    nAtom = size(iSquare)-1
    
    @:ASSERT(size(x) == size(y))
    @:ASSERT(mOrb >= 1)

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if
    
    if (present(beta)) then
      y = beta * y
    else
      y = 0.0_dp
    end if
    
    !$OMP PARALLEL DO DEFAULT(PRIVATE) REDUCTION(+:y) SCHEDULE(RUNTIME) &
    !$OMP & SHARED(nAtom, iSquare, nNeighbour, iNeighbour, img2CentCell, iPair)
    do iAtom1 = 1, nAtom
       nOrb1 = iSquare(iAtom1+1) - iSquare(iAtom1)
       ix = iSquare(iAtom1)
       jx = iSquare(iAtom1 + 1)-1
       @:ASSERT(nOrb1==jx-ix+1)
       do iNeigh = 0, nNeighbour(iAtom1)
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          nOrb2 = iSquare(iAtom2f+1)-iSquare(iAtom2f)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrTmp(1:nOrb2,1:nOrb1) = alphaTmp*reshape(A(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
          ! symmetrize on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii) 
            end do
          end if
          iy = iSquare(iAtom2f)
          jy = iSquare(iAtom2f + 1)-1
          @:ASSERT(nOrb2==jy-iy+1)
          y(iy:jy) = y(iy:jy) + matmul(sqrTmp(1:nOrb2,1:nOrb1),x(ix:jx))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            y(ix:jx) = y(ix:jx) + matmul(transpose(sqrTmp(1:nOrb2,1:nOrb1)),x(iy:jy))
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine SYMV_gamma


  !> Symmetric sparse matrix multiplied with a general matrix C = alpha * A * B + beta * C
  subroutine SYMM_gamma(C, A, B, iNeighbour, nNeighbour, img2CentCell, iPair, iSquare, mOrb, alpha,&
      & beta)

    !> Matrix to accumulate into
    real(dp), intent(inout) :: C(:,:)

    !> Symmetric matrix in sparse format
    real(dp), intent(in) :: A(:)

    !> General matrix
    real(dp), intent(in) :: B(:,:)

    !> Neighbouring atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbouring atoms
    integer, intent(in) :: nNeighbour(:)

    !> Atom indices into central cell atom numbers
    integer, intent(in) :: img2CentCell(:)

    !> Indexing for sparse matrix
    integer, intent(in) :: iPair(0:,:)

    !> Indexing for dense matrinx
    integer, intent(in) :: iSquare(:)

    !> Maximum number of orbitals in an atomic block
    integer, intent(in) :: mOrb

    !> Scaling factor for A B, 1 if not set
    real(dp), optional, intent(in) :: alpha

    !> Scaling factor for incomming C, 0 if not set
    real(dp), optional, intent(in) :: beta

    integer :: iOrig, iB, iC, jB, jC, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii
    real(dp) :: sqrTmp(mOrb,mOrb), alphaTmp
    
    nAtom = size(iSquare)-1
    
    @:ASSERT(all(shape(C) == shape(B)))
    @:ASSERT(mOrb >= 1)
    
    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if
    
    if (present(beta)) then
      C(:,:) = beta * C
    else
      C(:,:) = 0.0_dp
    end if
    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP& PRIVATE(iAtom1,nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC) &
    !$OMP& REDUCTION(+:C) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
       nOrb1 = iSquare(iAtom1+1)-iSquare(iAtom1)
       iB = iSquare(iAtom1)
       jB = iSquare(iAtom1 + 1)-1
       @:ASSERT(nOrb1==jB-iB+1)
       do iNeigh = 0, nNeighbour(iAtom1)
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          nOrb2 = iSquare(iAtom2f+1)-iSquare(iAtom2f)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrTmp(1:nOrb2,1:nOrb1) = alphaTmp*reshape(A(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
          ! symmetrize on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii) 
            end do
          end if
          iC = iSquare(iAtom2f)
          jC = iSquare(iAtom2f + 1)-1
          @:ASSERT(nOrb2==jC-iC+1)
          C(iC:jC, :) = C(iC:jC, :) + matmul(sqrTmp(1:nOrb2,1:nOrb1), B(iB:jB,:))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(iB:jB, :) = C(iB:jB, :) + matmul(transpose(sqrTmp(1:nOrb2,1:nOrb1)), B(iC:jC,:))
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine SYMM_gamma

end module dftbp_sparseblas
