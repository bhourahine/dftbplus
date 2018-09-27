!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#! Note: This module contains preprocessor variable substitutions in subroutine names (${NAME}$)
#! which may break the documentation system. Make sure you preprocess this file before passing it
#! to a source code documentation tool.

#:include 'common.fypp'

!> Interface wrapper for the blas routines.
!>
!> ALL BLAS routines which are called from the main code must be included here.
module blas
  use accuracy, only : rsp, rdp
  public

  interface

  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssyr'), ('rdp', 'real', 'dsyr'),&
    & ('rsp', 'complex', 'cher'), ('rdp', 'complex', 'zher')]
    !> performs the rank 1 operation
    !> A := alpha*x*x**T + A or A := alpha*x*x**H + A
    subroutine ${NAME}$(uplo, nn, alpha, xx, incx, aa, lda)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(${VPREC}$), intent(in) :: alpha

      !> vector
      ${VTYPE}$(${VPREC}$), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)
      
    end subroutine ${NAME}$
  #:endfor

    
  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sger'), ('rdp', 'real', 'dger'),&
    & ('rsp', 'complex', 'cgerc'), ('rdp', 'complex', 'zgerc')]
    !> performs the rank 1 operation
    !> A := alpha*x*y**T + A or A := alpha*x*y**H + A
    subroutine ${NAME}$(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import ${VPREC}$

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> vector
      ${VTYPE}$(${VPREC}$), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> vector
      ${VTYPE}$(${VPREC}$), intent(in) :: yy(*)

      !> stride
      integer, intent(in) :: incy

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)
      
    end subroutine ${NAME}$
  #:endfor

  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssymv'), ('rdp', 'real', 'dsymv'),&
    & ('rsp', 'complex', 'chemv'), ('rdp', 'complex', 'zhemv')]
    !> performs the matrix-vector operation
    !> y := alpha*A*x + beta*y
    subroutine ${NAME}$(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> vector
      ${VTYPE}$(${VPREC}$), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: beta

      !> vector
      ${VTYPE}$(${VPREC}$), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine ${NAME}$
  #:endfor
    

  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sgemv'), ('rdp', 'real', 'dgemv')]
    !> performs one of the matrix-vector operations
    !> y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
    subroutine ${NAME}$(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import ${VPREC}$

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> vector
      ${VTYPE}$(${VPREC}$), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: beta

      !> vector
      ${VTYPE}$(${VPREC}$), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
      
    end subroutine ${NAME}$
  #:endfor

    
  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssymm'), ('rdp', 'real', 'dsymm'),&
    & ('rsp', 'complex', 'chemm'), ('rdp', 'complex', 'zhemm')]
    !> performs one of the matrix-matrix operations
    !> C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
    subroutine ${NAME}$(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import ${VPREC}$

      !> side of the product
      character, intent(in) :: side

      !> upper or lower matrix
      character, intent(in) :: uplo

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(in) :: bb(ldb, *)

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      ${VTYPE}$(${VPREC}$), intent(inout) :: cc(ldc, *)
      
    end subroutine ${NAME}$
  #:endfor

    
  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sgemm'), ('rdp', 'real', 'dgemm'),&
    & ('rsp', 'complex', 'cgemm'), ('rdp', 'complex', 'zgemm')]
    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, op( X ) = X**T, or op( X ) = X**H
    subroutine ${NAME}$(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import ${VPREC}$

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(in) :: bb(ldb, *)

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      ${VTYPE}$(${VPREC}$), intent(inout) :: cc(ldc, *)
      
    end subroutine ${NAME}$
  #:endfor
    

  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssyrk'), ('rdp', 'real', 'dsyrk'),&
    & ('rsp', 'complex', 'cherk'), ('rdp', 'complex', 'zherk')]
    !> performs one of the symmetric rank k operations
    !> C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C
    subroutine ${NAME}$(uplo, trans, nn, kk, alpha, aa, lda, beta, cc, ldc)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix  size
      integer, intent(in) :: nn

      !> rank
      integer, intent(in) :: kk

      !> scaling factor
      real(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> scale factor
      real(${VPREC}$), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      ${VTYPE}$(${VPREC}$), intent(inout) :: cc(ldc, *)
      
    end subroutine ${NAME}$
  #:endfor

  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'strsm'), ('rdp', 'real', 'dtrsm'),&
    & ('rsp', 'complex', 'ctrsm'), ('rdp', 'complex', 'ztrsm')]
    !> Solves one of the matrix equations
    !> op( A )*X = alpha*B, or X*op( A ) = alpha*B
    subroutine ${NAME}$(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import ${VPREC}$

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(inout) :: bb(ldb, *)
      
    end subroutine ${NAME}$
  #:endfor

    
  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'strmm'), ('rdp', 'real', 'dtrmm'),&
    & ('rsp', 'complex', 'ctrmm'), ('rdp', 'complex', 'ztrmm')]
    !> Performs one of the matrix-matrix operations
    !> B := alpha*op( A )*B, or B := alpha*B*op( A ),
    !> where op( A ) = A, op( A ) = A**T or op( A ) = A**H.
    subroutine ${NAME}$(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import ${VPREC}$

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(inout) :: bb(ldb, *)
      
    end subroutine ${NAME}$
  #:endfor 

  end interface

end module blas
