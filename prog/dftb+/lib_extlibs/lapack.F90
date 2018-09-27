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

!> Interface wrapper for the lapack routines. See the <a href="http://www.netlib.org/lapack/">lapack
!> project documentation</a> for more details
module lapack
  use accuracy, only : rsp, rdp
  implicit none
  public

#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssyev'), ('rdp', 'real', 'dsyev'),&
  & ('rsp', 'complex', 'cheev'), ('rdp', 'complex', 'zheev')]
  !> Symmetric/hermitian eigensolver
  interface ${NAME}$

    !> Symmetric/hermitian eigensolver
    subroutine ${NAME}$(jobz, uplo, nn, aa, lda, ww, work, lwork,&
    #:if VTYPE == 'complex'
        & rwork, info)
    #:else
      & info)
    #:endif

      import ${VPREC}$

      !> job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(${VPREC}$), intent(out) :: ww(*)

      !> workspace
      ${VTYPE}$(${VPREC}$), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

    #:if VTYPE == 'complex'
      !> real workspace
      real(${VPREC}$), intent(inout) :: rwork(*)
    #:endif

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssygv'), ('rdp', 'real', 'dsygv'),&
  & ('rsp', 'complex', 'chegv'), ('rdp', 'complex', 'zhegv')]
  !> Simple generalised hermitian/symmetric eigensolver
  interface ${NAME}$

    !> Simple generalised hermitian/symmetric eigensolver
    subroutine ${NAME}$(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
      #:if VTYPE == 'complex'
        & rwork, info)
      #:else
        & info)
      #:endif

      import ${VPREC}$

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> leading dimension of B
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(${VPREC}$), intent(out) :: ww(*)

      !> workspace
      ${VTYPE}$(${VPREC}$), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

    #:if VTYPE == 'complex'
      !> real workspace
      real(${VPREC}$), intent(inout) :: rwork(*)
    #:endif

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssygvd'), ('rdp', 'real', 'dsygvd'),&
  & ('rsp', 'real', 'ssygvr'), ('rdp', 'real', 'dsygvr'),&
  & ('rsp', 'complex', 'chegvd'), ('rdp', 'complex', 'zhegvd'),&
  & ('rsp', 'complex', 'chegvr'), ('rdp', 'complex', 'zhegvr') ]
  !> Generalised hermitian/symmetric eigensolver
  interface ${NAME}$

    !> Generalised hermitian/symmetric eigensolver
    subroutine ${NAME}$(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
      #:if VTYPE == 'complex'
        & rwork, lrwork,&
      #:endif
        & iwork, liwork, info)

      import ${VPREC}$

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> leading dimension of B
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(${VPREC}$), intent(out) :: ww(*)

      !> workspace
      ${VTYPE}$(${VPREC}$), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

    #:if VTYPE == 'complex'
      !> real workspace
      real(${VPREC}$), intent(inout) :: rwork(*)

      !> size of rwork
      integer, intent(in) :: lrwork
    #:endif

      !> integer workspace
      integer, intent(inout) :: iwork(*)

      !> size of integer workspace
      integer, intent(in) :: liwork

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'spotrf'), ('rdp', 'real', 'dpotrf'),&
  & ('rsp', 'complex', 'cpotrf'), ('rdp', 'complex', 'zpotrf')]
  !> Cholesky factorization of symmetric/hermitian positive definite matrix
  interface ${NAME}$

    !> Cholesky factorization of double ${VTYPE}$ hermitian positive definite matrix
    subroutine ${NAME}$(uplo, nn, aa, lda, info)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssygst'), ('rdp', 'real', 'dsygst'),&
  & ('rsp', 'complex', 'chegst'), ('rdp', 'complex', 'zhegst')]
  !> Reduce a symmetric/hermitian definite generalized eigenproblem to standard form
  interface ${NAME}$

    !> Reduce double complex hermitian-definite generalized eigenproblem to standard form
    subroutine ${NAME}$(itype, uplo, nn, aa, lda, bb, ldb, info)
      import ${VPREC}$

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> leading dimension of B
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(in) :: bb(ldb, *)

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sgesv'), ('rdp', 'real', 'dgesv')]
  !> Solve overdetermined or underdetermined linear systems
  interface ${NAME}$

    !> Solve overdetermined or underdetermined real linear systems
    subroutine ${NAME}$(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import ${VPREC}$

      !> matrix dimension
      integer, intent(in) :: nn

      !> number of right hand side equations
      integer, intent(in) :: nrhs

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> pivot array
      integer, intent(out) :: ipiv(*)

      !> leading dimension of B
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(inout) :: bb(ldb, *)

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sgetrf'), ('rdp', 'real', 'dgetrf'),&
  & ('rsp', 'complex', 'cgetrf'), ('rdp', 'complex', 'zgetrf')]
  !> Computes LU factorization of general matrix
  interface ${NAME}$

    !> Computes LU factorization of general matrix
    subroutine ${NAME}$(mm, nn, aa, lda, ipiv, info)
      import ${VPREC}$

      !> number of rows of the matrix
      integer, intent(in) :: mm

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> pivot array
      integer, intent(out) :: ipiv(*)

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sgetri'), ('rdp', 'real', 'dgetri')]
  !> Computes inverse of a general matrix using LU factorisation
  interface ${NAME}$

    !> Computes inverse of a general matrix using LU factorisation
    subroutine ${NAME}$(nn, aa, lda, ipiv, work, lwork, info)
      import ${VPREC}$

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> pivot array
      integer, intent(in) :: ipiv(*)

      !> workspace
      ${VTYPE}$(${VPREC}$), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssytrf'), ('rdp', 'real', 'dsytrf')]
  !> Factorise a symmetric matrix as A = U*D*U**T or A = L*D*L**T
  interface ${NAME}$

    !> Factorise a symmetric matrix as A = U*D*U**T or A = L*D*L**T
    subroutine ${NAME}$(uplo, nn, aa, lda, ipiv, work, lwork, info)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(inout) :: aa(lda, *)

      !> pivot array
      integer, intent(out) :: ipiv(*)

      !> workspace
      ${VTYPE}$(${VPREC}$), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssytrs'), ('rdp', 'real', 'dsytrs')]
  !> Solve a system of linear equations for a symmetric matrix A*X = B
  interface ${NAME}$

    !> Solve a system of linear equations for a symmetric matrix A*X = B
    subroutine ${NAME}$(uplo, nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> number of right hand side equations
      integer, intent(in) :: nrhs

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> pivot array
      integer, intent(in) :: ipiv(*)

      !> leading dimension of B
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(inout) :: bb(ldb, *)

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'slarnv'), ('rdp', 'real', 'dlarnv'),&
  & ('rsp', 'complex', 'clarnv'), ('rdp', 'complex', 'zlarnv')]
  !> Returns a vector of random numbers from a uniform or normal distribution
  interface ${NAME}$

    !> Returns a vector of random numbers from a uniform or normal distribution
    subroutine ${NAME}$(idist, iseed, nn, xx)
      import ${VPREC}$

      !> distribution choice
      integer, intent(in) :: idist

      !> generator seed
      integer, intent(inout) :: iseed(4)

      !> vector dimension
      integer, intent(in) :: nn

      !> Random values on exit
      ${VTYPE}$(${VPREC}$), intent(out) :: xx(*)

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor


  !> Provides problem-dependent LAPACK routine parameters for the local environment
  interface ilaenv

    !> Provides problem-dependent LAPACK routine parameters for the local environment
    function ilaenv(ispec, name, opts, n1, n2, n3, n4)

      !> Specifies the parameter to be returned
      integer, intent(in) :: ispec

      !> name of alling subroutine
      character, intent(in) :: name

      !> The character options to the subroutine NAME, concatenated together
      character, intent(in) :: opts

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n1

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n2

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n3

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n4

      !> returned parameter
      integer :: ilaenv

    end function ilaenv

  end interface ilaenv


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'slamch'), ('rdp', 'real', 'dlamch')]
  !> Machine parameters
  interface ${NAME}$

    !> Machine parameters
    function ${NAME}$(cmach)
      import ${VPREC}$

      !> name of parameter to return
      character, intent(in) :: cmach

      !> parameter value
      ${VTYPE}$(${VPREC}$) :: ${NAME}$
    end function ${NAME}$

  end interface ${NAME}$
#:endfor


  !> Error handler for the LAPACK routines
  interface xerbla

    !> Error handler for the LAPACK routines
    subroutine xerbla(srname, info)

      !> calling subroutine name
      character(6), intent(in) :: srname

      !> info state of the routine
      integer, intent(in) :: info
    end subroutine xerbla
  end interface xerbla


#:for VPREC, VTYPE, NAME in [('rsp', 'real', 'sgesvd'), ('rdp', 'real', 'dgesvd'),&
  & ('rsp', 'complex', 'cgesvd'), ('rdp', 'complex', 'zgesvd')]
  !> Singular value decomposition
  interface ${NAME}$

    !> Singular value decomposition
    subroutine ${NAME}$(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork,&
    #:if VTYPE == 'complex'
        & rwork, info)
    #:else
        & info)
    #:endif

      import ${VPREC}$

      !> job type for vt
      character, intent(in) :: jobvt

      !> job type for u
      character, intent(in) :: jobu

      !> First matrix dimension for A
      integer, intent(in) :: m

      !> Second matrix dimension for A
      integer, intent(in) :: n

      !> leading dimension of A
      integer, intent(in) :: lda

      !> leading dimension of U
      integer, intent(in) :: ldu

      !> leading dimension of Vt
      integer, intent(in) :: ldvt

      !> matrix to decompose
      ${VTYPE}$(${VPREC}$), intent(inout) :: a(lda,*)

      !> singular values on return min(m,n)
      real(${VPREC}$), intent(out) :: s(*)

    #:if VTYPE == 'complex'
      !> real workspace
      real(${VPREC}$), intent(out) :: rwork(*)
    #:endif

      !> Left singular vectors
      ${VTYPE}$(${VPREC}$), intent(out) :: u(ldu,*)

      !> Right singular vectors
      ${VTYPE}$(${VPREC}$), intent(out) :: vt(ldvt,*)

      !> ${VTYPE}$ work space
      ${VTYPE}$(${VPREC}$), intent(out) :: work(*)

      !> size of ${VTYPE}$ work space
      integer, intent(in) :: lwork

      !> state of routine on return
      integer, intent(in) :: info

    end subroutine ${NAME}$

  end interface ${NAME}$
#:endfor

end module lapack
