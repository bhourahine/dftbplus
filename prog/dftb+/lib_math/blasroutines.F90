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

!> Contains F90 wrapper functions for some commonly used blas calls needed in the code. The
!> interface of all BLAS calls must be defined in the module blas.
module blasroutines
  use assert
  use accuracy
  use blas
  implicit none


  !> Rank 1 update of a matrix A := alpha*x*x' + A
  !> Wrapper for the level 2 blas routine xsyr to perform the rank 1 update of the chosen triangle
  interface her
  #:for IFACETYPE in [('real'), ('cmplx'), ('dble'), ('dblecmplx')]
    module procedure her_${IFACETYPE}$
  #:endfor
  end interface her


  !> Rank 1 update of a matrix A := alpha*x*y' + A
  !> Wrapper for the level 2 blas routine xger to perform the rank 1 update of the chosen triangle
  interface ger
  #:for IFACETYPE in [('real'), ('cmplx'), ('dble'), ('dblecmplx')]
    module procedure ger_${IFACETYPE}$
  #:endfor
  end interface ger


  !> Symmetric matrix vector multiply y := alpha*A*x + beta*y
  !> Wrapper for the level 2 blas routine
  interface hemv
  #:for IFACETYPE, NAME in [('real', 'symv'), ('cmplx', 'hemv'), ('dble', 'symv'),&
    & ('dblecmplx', 'hemv')]
    module procedure ${NAME}$_${IFACETYPE}$
  #:endfor
  end interface hemv


  !> General matrix vector multiply y := alpha*a*x + beta*y
  !> Wrapper for the level 2 blas routine
  interface gemv
    module procedure gemv_real
    module procedure gemv_dble
  end interface gemv


  !> Interface to GEMM routines evaluates C := alpha*op( A )*op( B ) + beta*C, where op( X ) is one
  !> of op( X ) = X or op( X ) = X'
  interface gemm
  #:for IFACETYPE in [('real'), ('cmplx'), ('dble'), ('dblecmplx')]
    module procedure gemm_${IFACETYPE}$
  #:endfor
  end interface gemm


  !> Rank k update of a matrix C := alpha*A*A' + beta C
  !> Wrapper for the level 3 blas routine syrk/herk to perform the rank k update of the chosen
  !> triangle of C
  interface herk
  #:for IFACETYPE in [('real'), ('cmplx'), ('dble'), ('dblecmplx')]
    module procedure herk_${IFACETYPE}$
  #:endfor
  end interface herk


  !> Interface to SYMM routines
  !> Wrapper for the level 3 blas routine
  interface symm
    module procedure symm_real
    module procedure symm_dble
  end interface symm


  !> Interface to HEMM routines
  !> Wrapper for the level 3 blas routine
  interface hemm
    module procedure hemm_cmplx
    module procedure hemm_dblecmplx
  end interface hemm


contains


#:for LABEL, VTYPE, VPREC, NAME in [('real', 'real', 'rsp', 'ssyr'),&
  & ('cmplx', 'complex', 'rsp', 'cher'), ('dble', 'real', 'rdp', 'dsyr'),&
  & ('dblecmplx', 'complex', 'rdp', 'zher')]
  !> Rank 1 update of a symmetric/hermitian matrix
  subroutine her_${LABEL}$(a,alpha,x,uplo)

    !> contains the matrix for the update
    ${VTYPE}$(${VPREC}$), intent(inout) :: a(:,:)

    !> scaling value for the update contribution
    real(${VPREC}$), intent(in) :: alpha

    !> vector of values for the update
    ${VTYPE}$(${VPREC}$), intent(in) :: x(:)

    !> optional upper, 'U', or lower 'L' triangle, defaults to lower
    character, intent(in),optional :: uplo

    integer :: n
    character :: iuplo

    @:ASSERT(size(a,dim=1) == size(a,dim=2))
    @:ASSERT(size(a,dim=1) == size(x))

    if (present(uplo)) then
      iuplo = uplo
    else
      iuplo = 'l'
    end if
    @:ASSERT(iuplo == 'u' .or. iuplo == 'U' .or. iuplo == 'l' .or. iuplo == 'L')

    n = size(x)
    call ${NAME}$(iuplo,n,alpha,x,1,a,n)

  end subroutine her_${LABEL}$
#:endfor


#:for LABEL, VTYPE, VPREC, NAME in [('real', 'real', 'rsp', 'sger'),&
  & ('cmplx', 'complex', 'rsp', 'cgerc'), ('dble', 'real', 'rdp', 'dger'),&
  & ('dblecmplx', 'complex', 'rdp', 'zgerc')]
  !> Rank 1 update of a general matrix
  subroutine ger_${LABEL}$(a, alpha, x, y)

    !> contains the matrix for the update
    ${VTYPE}$(${VPREC}$), intent(inout) :: a(:,:)

    !> scaling value for the update contribution
    ${VTYPE}$(${VPREC}$), intent(in) :: alpha

    !> vector of values for the update
    ${VTYPE}$(${VPREC}$), intent(in) :: x(:)

    !> vector of values for the update
    ${VTYPE}$(${VPREC}$), intent(in) :: y(:)

    integer :: n, m

    @:ASSERT(size(a,dim=1) == size(x))
    @:ASSERT(size(a,dim=2) == size(y))

    m = size(x)
    n = size(y)

    call ${NAME}$(m, n, alpha, x, 1, y, 1, a, m)

  end subroutine ger_${LABEL}$
#:endfor


#:for LABEL, VTYPE, VPREC, NAME, CASE in [('real', 'real', 'rsp', 'symv', 's'),&
  & ('cmplx', 'complex', 'rsp', 'hemv', 'c'), ('dble', 'real', 'rdp', 'symv', 'd'),&
  & ('dblecmplx', 'complex', 'rdp', 'hemv', 'z')]
  !> Symmetric/hermitian matrix*vector product
  subroutine ${NAME}$_${LABEL}$(y, a, x, uplo, alpha, beta)

    !> vector
    ${VTYPE}$(${VPREC}$), intent(inout) :: y(:)

    !> symmetric matrix
    ${VTYPE}$(${VPREC}$), intent(in) :: a(:,:)

    !> vector
    ${VTYPE}$(${VPREC}$), intent(in) :: x(:)

    !> optional upper, 'U', or lower 'L' triangle (defaults to lower)
    character, intent(in), optional :: uplo

    !> optional scaling factor (defaults to 1)
    ${VTYPE}$(${VPREC}$), intent(in), optional :: alpha

    !> optional scaling factor (defaults to 0)
    ${VTYPE}$(${VPREC}$), intent(in), optional :: beta

    integer :: n
    character :: iUplo
    ${VTYPE}$(${VPREC}$) :: iAlpha, iBeta

    if (present(uplo)) then
      iUplo = uplo
    else
      iUplo = 'L'
    end if
    if (present(alpha)) then
      iAlpha = alpha
    else
    #:if VTYPE == 'complex'
      iAlpha = cmplx(1.0,0.0,${VPREC}$)
    #:else
      iAlpha = 1.0_${VPREC}$
    #:endif
    end if
    if (present(beta)) then
      iBeta = beta
    else
    #:if VTYPE == 'complex'
      iBeta = cmplx(0.0,0.0,${VPREC}$)
    #:else
      iBeta = 0.0_${VPREC}$
    #:endif
    end if

    @:ASSERT(size(y) == size(x))
    @:ASSERT(size(a,dim=1) == size(a,dim=2))
    @:ASSERT(size(a,dim=1) == size(x))
    @:ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')

    n = size(y)

    call ${CASE}$${NAME}$(iUplo, n, iAlpha, a, n, x, 1, iBeta, y, 1)

  end subroutine ${NAME}$_${LABEL}$
#:endfor


#:for IFACETYPE, VPREC, VTYPE in [('real', 'rsp', 's'), ('dble', 'rdp', 'd')]
  !> General matrix*vector produc
  subroutine gemv_${IFACETYPE}$(y, a, x, alpha, beta, trans)

    !> vector
    real(${VPREC}$), intent(inout) :: y(:)

    !> matrix
    real(${VPREC}$), intent(in) :: a(:,:)

    !> vector
    real(${VPREC}$), intent(in) :: x(:)

    !> optional scaling factor (defaults to 1)
    real(${VPREC}$), intent(in), optional :: alpha

    !> optional scaling factor (defaults to 0)
    real(${VPREC}$), intent(in), optional :: beta

    !> optional transpose (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c' and 'C'
    character, intent(in), optional :: trans

    integer :: n, m
    character :: iTrans
    real(${VPREC}$) :: iAlpha, iBeta

    if (present(trans)) then
      iTrans = trans
    else
      iTrans = 'n'
    end if
    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = 1.0_${VPREC}$
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = 0.0_${VPREC}$
    end if

    @:ASSERT(iTrans == 'n' .or. iTrans == 'N' .or. iTrans == 't' .or. iTrans == 'T' .or.&
        & iTrans == 'c' .or. iTrans == 'C')
    @:ASSERT(((size(a,dim=1) == size(y)) .and. (iTrans == 'n' .or. iTrans == 'N')) .or.&
        & (size(a,dim=1) == size(x)))
    @:ASSERT(((size(a,dim=2) == size(x)) .and. (iTrans == 'n' .or. iTrans == 'N')) .or.&
        & (size(a,dim=2) == size(y)))

    m = size(a,dim=1)
    n = size(a,dim=2)

    call ${VTYPE}$gemv( iTrans, m, n, iAlpha, a, m, x, 1, iBeta, y, 1 )

  end subroutine gemv_${IFACETYPE}$
#:endfor


#:for IFACETYPE, VPREC, VTYPE in [('real', 'rsp', 's'), ('dble', 'rdp', 'd')]
  !> Real symmetric matrix * general matrix multiply
  subroutine symm_${IFACETYPE}$(C, side, A, B, uplo, alpha, beta, m, n)

    !> general matrix output
    real(${VPREC}$), intent(inout) :: C(:,:)

    !> symmetric matrix on 'l'eft or 'r'ight , where A is symmetric and B is general SIDE = 'L' or
    !> 'l' C := alpha*A*B + beta*C, SIDE = 'R' or 'r' C := alpha*B*A + beta*C
    character, intent(in) :: side

    !> symmetric matrix, size
    real(${VPREC}$), intent(in) :: A(:,:)

    !> general matrix
    real(${VPREC}$), intent(in) :: B(:,:)

    !> is an 'U'pper or 'L'ower triangle matrix, defaults to lower
    character, intent(in), optional :: uplo

    !> defaults to 1 if not set
    real(${VPREC}$), intent(in), optional :: alpha

    !> defaults to 0 if not set
    real(${VPREC}$), intent(in), optional :: beta

    !> specifies the number of rows of the matrix C
    integer, intent(in), optional :: m

    !> specifies the number of columns of the matrix C
    integer, intent(in), optional :: n

    integer :: lda, ldb, ldc, ka, im, in
    character :: iUplo
    real(${VPREC}$) :: iAlpha, iBeta

    if (present(uplo)) then
      iUplo = uplo
    else
      iUplo = 'L'
    end if
    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = 1.0_${VPREC}$
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = 0.0_${VPREC}$
    end if

    lda = size(a,dim=1)
    ldb = size(b,dim=1)
    ldc = size(c,dim=1)

    @:ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')
    @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')

    if (present(n)) then
      in = n
    else
      in = size(C,dim=2)
    end if
    if (present(m)) then
      im = m
    else
      im = size(C,dim=1)
    end if
    if (iUplo=='l' .or. iUplo=='L') then
      ka = im
    else
      ka = in
    end if

    @:ASSERT(im>0)
    @:ASSERT(in>0)
    @:ASSERT(lda>=im .and. (iUplo=='l' .or. iUplo=='L') .or. lda>=in)
    @:ASSERT(size(a,dim=2)>=ka)
    @:ASSERT(ldb>=im)
    @:ASSERT(ldc>=im)
    @:ASSERT(size(B,dim=2)>=in)
    @:ASSERT(size(C,dim=2)>=in)

    call ${VTYPE}$symm(side, iUplo, im, in, iAlpha, A, lda, B, ldb, iBeta, C, ldc)

  end subroutine symm_${IFACETYPE}$
#:endfor


#:for LABEL, VTYPE, VPREC, CASE in [('real', 'real', 'rsp', 's'),&
  & ('cmplx', 'complex', 'rsp', 'c'), ('dble', 'real', 'rdp', 'd'),&
  & ('dblecmplx', 'complex', 'rdp', 'z')]
  !> Matrix*matrix product
  subroutine gemm_${LABEL}$(C, A, B, alpha, beta, transA, transB, n, m, k)

    !> general matrix output
    ${VTYPE}$(${VPREC}$), intent(inout) :: C(:,:)

    !> symmetric matrix
    ${VTYPE}$(${VPREC}$), intent(in) :: A(:,:)

    !> general matrix
    ${VTYPE}$(${VPREC}$), intent(in) :: B(:,:)

    !> defaults to 1 if not set
    ${VTYPE}$(${VPREC}$), intent(in), optional :: alpha

    !> defaults to 0 if not set
    ${VTYPE}$(${VPREC}$), intent(in), optional :: beta

    !> optional transpose of A matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c'
    !> and 'C'
    character, intent(in), optional :: transA

    !> optional transpose of B matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c'
    !> and 'C'
    character, intent(in), optional :: transB

    !> specifies the number of columns of the matrix C
    integer, intent(in), optional :: n

    !> specifies the number of rows of the matrix C
    integer, intent(in), optional :: m

    !> specifies the internal number of elements in Op(A)_ik Op(B)_kj
    integer, intent(in), optional :: k

    integer :: lda, ldb, ldc
    integer :: in, im, ik
    character :: iTransA, iTransB
    ${VTYPE}$(${VPREC}$) :: iAlpha, iBeta

    if (present(transA)) then
      iTransA = transA
    else
      iTransA = 'n'
    end if
    if (present(transB)) then
      iTransB = transB
    else
      iTransB = 'n'
    end if

    @:ASSERT(iTransA == 'n' .or. iTransA == 'N' .or. iTransA == 't'&
        & .or. iTransA == 'T' .or. iTransA == 'c' .or. iTransA == 'C')
    @:ASSERT(iTransB == 'n' .or. iTransB == 'N' .or. iTransB == 't'&
        & .or. iTransB == 'T' .or. iTransB == 'c' .or. iTransB == 'C')

    if (present(alpha)) then
      iAlpha = alpha
    else
    #:if VTYPE == 'complex'
      iAlpha = cmplx(1.0,0.0,${VPREC}$)
    #:else
      iAlpha = 1.0_${VPREC}$
    #:endif
    end if
    if (present(beta)) then
      iBeta = beta
    else
    #:if VTYPE == 'complex'
      iBeta = cmplx(0.0,0.0,${VPREC}$)
    #:else
      iBeta = 0.0_${VPREC}$
    #:endif
    end if

    lda = size(a,dim=1)
    ldb = size(b,dim=1)
    ldc = size(c,dim=1)

    if (present(m)) then
      im = m
    else
      if (iTransA == 'n' .or. iTransA == 'N') then
        im = size(A,dim=1)
      else
        im = size(A,dim=2)
      end if
    end if
    if (present(n)) then
      in = n
    else
      in = size(c,dim=2)
    end if
    if (present(k)) then
      ik = k
    else
      if (iTransA == 'n' .or. iTransA == 'N') then
        ik = size(A,dim=2)
      else
        ik = size(A,dim=1)
      end if
    end if

    @:ASSERT(im>0)
    @:ASSERT(in>0)
    @:ASSERT(ik>0)
    @:ASSERT((lda>=im .and. (iTransA == 'n' .or. iTransA == 'N')) .or. size(a,dim=2)>=im)
    @:ASSERT(ldc>=im)
    @:ASSERT((size(b,dim=2)>=in .and. (iTransB == 'n' .or. iTransB == 'N')) .or. ldb>=in)
    @:ASSERT(size(c,dim=2)>=in)
    @:ASSERT((size(a,dim=2)>=ik .and. (iTransA == 'n' .or. iTransA == 'N')) .or. lda>=ik)
    @:ASSERT((ldb>=ik .and. (iTransB == 'n' .or. iTransB == 'N')) .or. size(b,dim=2)>=ik)

    call ${CASE}$gemm(iTransA, iTransB, im, in, ik, iAlpha, A, lda, B, ldb, iBeta, C, ldc)

  end subroutine gemm_${LABEL}$
#:endfor


#:for LABEL, VTYPE, VPREC, NAME in [('real', 'real', 'rsp', 'ssyrk'),&
  & ('cmplx', 'complex', 'rsp', 'cherk'), ('dble', 'real', 'rdp', 'dsyrk'),&
  & ('dblecmplx', 'complex', 'rdp', 'zherk')]
  !> Rank-k update
  subroutine herk_${LABEL}$(C,A,alpha,beta,uplo,trans,n,k)

    !> contains the matrix to be updated
    ${VTYPE}$(${VPREC}$), intent(inout) :: C(:,:)

    !> contains the matrix to update
    ${VTYPE}$(${VPREC}$), intent(in) :: A(:,:)

    !> scaling value for the update contribution, defaults to 1
    real(${VPREC}$), intent(in), optional :: alpha

    !> scaling value for the original C, defaults to 0
    real(${VPREC}$), intent(in), optional :: beta

    !> optional upper, 'U', or lower 'L' triangle, defaults to lower
    character, intent(in), optional :: uplo

    !> optional transpose (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T' (and 'C' or 'c'
    !> for the real cases)
    character, intent(in), optional :: trans

    !> order of the matrix C
    integer, intent(in), optional :: n

    !> internal order of A summation
    integer, intent(in), optional :: k

    integer :: lda, ldc
    integer :: in, ik
    character :: iTrans, iUplo
    real(${VPREC}$) :: iAlpha, iBeta

    if (present(uplo)) then
      iUplo = uplo
    else
      iUplo = 'L'
    end if
    @:ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')

    if (present(trans)) then
      iTrans = trans
    else
      iTrans = 'n'
    end if

    @:ASSERT(iTrans == 'n' .or. iTrans == 'N' .or. iTrans == 't' .or. iTrans == 'T')

    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = 1.0_${VPREC}$
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = 0.0_${VPREC}$
    end if

    lda = size(a,dim=1)
    ldc = size(c,dim=1)

    if (present(n)) then
      in = n
    else
      in = size(c,dim=2)
    end if
    if (present(k)) then
      ik = k
    else
      if (iTrans == 'n' .or. iTrans == 'N') then
        ik = size(A,dim=2)
      else
        ik = size(A,dim=1)
      end if
    end if

    @:ASSERT(in>0)
    @:ASSERT(ik>0)
    @:ASSERT(((size(a,dim=2)>=in).and.(iTrans == 'n' .or. iTrans == 'N')) .or. (lda>=in))
    @:ASSERT(size(c,dim=2)>=in)
    @:ASSERT(((size(a,dim=2)>=ik).and.(iTrans == 'n' .or. iTrans == 'N')) .or. (lda>=ik))

    call ${NAME}$(iUplo, iTrans, in, ik, iAlpha, A, lda, iBeta, C, ldc )

  end subroutine herk_${LABEL}$
#:endfor


#:for LABEL, VPREC, CASE in [('cmplx', 'rsp', 'c'), ('dblecmplx', 'rdp', 'z')]
  !> double precision hermitian matrix * general matrix multiply
  subroutine hemm_${LABEL}$(C, side, A, B, uplo, alpha, beta, m, n)

    !> general matrix output
    complex(${VPREC}$), intent(inout) :: C(:,:)

    !> symmetric matrix on 'l'eft or 'r'ight , where A is symmetric and B is gen
    !> 'l' C := alpha*A*B + beta*C, SIDE = 'R' or 'r' C := alpha*B*A + beta*C
    character, intent(in) :: side

    !> hermitian matrix
    complex(${VPREC}$), intent(in) :: A(:,:)

    !> general matrix
    complex(${VPREC}$), intent(in) :: B(:,:)

    !> A is an 'U'pper or 'L'ower triangle matrix, defaults to lower
    character, intent(in), optional :: uplo

    !> defaults to 1 if not set
    complex(${VPREC}$), intent(in), optional :: alpha

    !> defaults to 0 if not set
    complex(${VPREC}$), intent(in), optional :: beta

    !> specifies the number of rows of the matrix C
    integer, intent(in), optional :: m

    !> specifies the number of columns of the matrix C
    integer, intent(in), optional :: n

    integer :: lda, ldb, ldc, ka, im, in
    character :: iUplo
    complex(${VPREC}$) :: iAlpha, iBeta

    if (present(uplo)) then
      iUplo = uplo
    else
      iUplo = 'L'
    end if
    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = (1.0_${VPREC}$, 0.0_${VPREC}$)
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = (0.0_${VPREC}$, 0.0_${VPREC}$)
    end if

    lda = size(a,dim=1)
    ldb = size(b,dim=1)
    ldc = size(c,dim=1)

    @:ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')
    @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')

    if (present(n)) then
      in = n
    else
      in = size(C,dim=2)
    end if
    if (present(m)) then
      im = m
    else
      im = size(C,dim=1)
    end if
    if (iUplo=='l'.or.iUplo=='L') then
      ka = im
    else
      ka = in
    end if

    @:ASSERT(im>0)
    @:ASSERT(in>0)
    @:ASSERT((lda>=im).and.(iUplo=='l'.or.iUplo=='L').or.(lda>=in))
    @:ASSERT(size(a,dim=2)>=ka)
    @:ASSERT(ldb>=im)
    @:ASSERT(ldc>=im)
    @:ASSERT(size(B,dim=2)>=in)
    @:ASSERT(size(C,dim=2)>=in)

    call ${CASE}$hemm(side, iUplo, im, in, iAlpha, A, lda, B, ldb, iBeta, C, ldc)

  end subroutine hemm_${LABEL}$
#:endfor

end module blasroutines
