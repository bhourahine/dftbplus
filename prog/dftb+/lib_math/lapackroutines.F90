!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains F90 wrapper functions for some commonly used lapack calls needed in the code. The
!> interface of all LAPACK calls must be defined in the module lapack.
module lapackroutines
  use assert
  use accuracy
  use message
  use lapack
  implicit none

  private


  !> Used to return runtime diagnostics
  character(len=100) :: error_string


  !> Computes the solution to a real system of linear equations A * X = B, where A is an N-by-N
  !> matrix and X and B are N-by-NRHS matrices
  interface gesv
    module procedure gesv_real
    module procedure gesv_dble
  end interface gesv


  !> Computes the LU decomposition of a general rectangular matrix using partial pivoting with row
  !> interchanges.
  !> The decomposition has the form: A = P*L*U, where P is a permutation matrix, L is a lower
  !> triangular matrix with unit diagonal elements and U is an upper triangular matrix.
  interface getrf
#:for NAME in [('real'), ('dble'), ('complex'), ('dcomplex')]
    module procedure getrf_${NAME}$
#:endfor
  end interface getrf

#:for NAME, VPREC in [('sytrf', 'real'), ('sytrf', 'dreal'), ('hetrf', 'complex'),&
  & ('hetrf', 'dcomplex')]
#:if NAME == 'sytrf'
  !> Bunch-Kaufman factorization of a symmetric matrix.
#:else
  !> Bunch-Kaufman factorization of a hermitian matrix.
#:endif
  interface ${NAME}$
    module procedure ${NAME}$_${VPREC}$
  end interface ${NAME}$
#:endfor


#:for NAME, VPREC in [('sytri', 'real'), ('sytri', 'dreal'), ('hetri', 'complex'),&
  & ('hetri', 'dcomplex')]
#:if NAME == 'sytri'
  !> Inverts a symmetric matrix.
#:else
  !> Inverts a Hermitian matrix.
#:endif
  interface ${NAME}$
    module procedure ${NAME}$_${VPREC}$
  end interface ${NAME}$
#:endfor


  !> Computes the inverse of a matrix using LU factorization computed by getrf.
  interface getri
    module procedure getri_real
    module procedure getri_dble
  end interface getri


  !> Solves a system of linear equations A*X = B with a real symmetric matrix A using the
  !> factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF.
  interface sytrs
    module procedure sytrs_dble
    module procedure sytrs_real
  end interface sytrs


  !> returns a vector of random numers, either from a uniform or normal distribution
  interface larnv
#:for NAME in [('real'), ('dble'), ('cplx'), ('dblecplx')]
    module procedure larnv_${NAME}$
#:endfor
  end interface larnv

  !> svd decomposition of matrix A into left and right vectors and singular values U S V^dag
  interface gesvd
#:for NAME, VPREC in [('sgesvd', 'real'), ('dgesvd', 'dble'), ('cgesvd', 'cplx'),&
    & ('zgesvd', 'dblecplx')]
    module procedure ${NAME}$_${VPREC}$
#:endfor
  end interface gesvd


  public :: gesv, getri, getrf, sytri, sytrf, matinv, symmatinv, sytrs, larnv
  public :: hermatinv, hetri, hetrf, gesvd


contains

#:for NAME, VPREC, VTYPE in [('real', 'rsp', 's'), ('dble', 'rdp', 'd')]
  !> Computes the solution to a real system of linear equations A * X = B
  subroutine gesv_${NAME}$(aa, bb, nEquation, nSolution, iError)

    !> Contains the coefficients on entry, the LU factorisation on exit.
    real(${VPREC}$), intent(inout) :: aa(:,:)

    !> Right hand side(s) of the linear equation on entry, solution(s) on exit.
    real(${VPREC}$), intent(inout) :: bb(:,:)

    !> The size of the problem (nr. of variables and equations). Must be only specified if different
    !> from size(aa, dim=1).
    integer, intent(in), optional :: nEquation

    !> Nr. of right hand sides (nr. of solutions). Must be only specified if different from size(b,
    !> dim=2).
    integer, intent(in), optional :: nSolution

    !> Error flag. If present, Lapack error flags are reported and noncritical errors (iError > 0)
    !> will not abort the program.
    integer, intent(out), optional :: iError
    integer :: info
    integer :: nn, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)

    lda = size(aa, dim=1)
    if (present(nEquation)) then
      @:ASSERT(nEquation >= 1 .and. nEquation <= lda)
      nn = nEquation
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)

    ldb = size(bb, dim=1)
    @:ASSERT(ldb >= nn)
    nrhs = size(bb, dim=2)
    if (present(nSolution)) then
      @:ASSERT(nSolution <= nrhs)
      nrhs = nSolution
    end if

    info = 0
    allocate(ipiv(nn))
    call ${VTYPE}$gesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)

    if (info < 0) then
      write (error_string, "('Failure in linear equation solver ${VTYPE}$gesv, illegal argument at&
          & position ',I2)") info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
        write (error_string, "('Linear dependent system in linear equation solver ${VTYPE}$getrf,&
            & info flag: ',I0)") info
        call error(error_string)
      end if
    end if

  end subroutine gesv_${NAME}$
#:endfor


#:for NAME, VPREC, VTYPE, CASE in [('real', 'rsp', 'real', 's'), ('dble', 'rdp', 'real', 'd'),&
  & ('complex', 'rsp', 'complex', 'c'), ('dcomplex', 'rdp', 'complex', 'z')]
  !> Computes the LU decomposition of a general rectangular matrix
  subroutine getrf_${NAME}$(aa, ipiv, nRow, nColumn, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    ${VTYPE}$(${VPREC}$), intent(inout) :: aa(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(out) :: ipiv(:)

    !> Number of rows of the matrix to decomposea. (Necessary if different from the number of rows
    !> of the passed matrix)
    integer, optional, intent(in) :: nRow

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of columns
    !> of the passed matrix)
    integer, optional, intent(in) :: nColumn

    !> Error flag. Zero on successfull exit. If not present, any lapack error causes program
    !> termination. If passed only fatal lapack errors with error flag < 0 cause abort.
    integer, optional, intent(out) :: iError

    integer :: mm, nn, lda, info

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      @:ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    @:ASSERT(size(ipiv) == min(mm, nn))

    call ${CASE}$getrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
      write (error_string, "('Failure in LU factorisation ${CASE}$getrf, illegal argument at&
          & position ',I2)") info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
        write (error_string, "('Factor U is exactly zero in ${CASE}$getrf, info flag is ',i10)")&
            & info
        call error(error_string)
      end if
    end if

  end subroutine getrf_${NAME}$
#:endfor


#:for NAME, VPREC, VTYPE in [('real', 'rsp', 's'), ('dble', 'rdp', 'd')]
  !> Computes the inverse of a matrix using LU factorization
  subroutine getri_${NAME}$(aa, ipiv, nRow, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    real(${VPREC}$), intent(inout) :: aa(:,:)

    !> Pivot indices, as calculated by getri
    integer, intent(in) :: ipiv(:)

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of rows of
    !> the passed matrix)
    integer, intent(in), optional :: nRow

    !> iError Error flag. Zero on successfull exit. If not present, any lapack error causes program
    !> termination. If present, only fatal lapack errors with error flag < 0 cause abort.
    integer, intent(out), optional :: iError

    integer :: nn, lda, info, lwork
    real(${VPREC}$), allocatable :: work(:)
    real(${VPREC}$) :: work2(1)

    lda = size(aa, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)
    @:ASSERT(size(ipiv) == nn)

    ! workspace query
    lwork = -1
    call ${VTYPE}$getri(nn, aa, lda, ipiv, work2, lwork, info)
    if (info < 0) then
      write (error_string, "('Failure in LU factorisation memory query (${VTYPE}$getri), illegal&
          & argument at position ',I2)") info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
        write (error_string, "('Failure in LU factorisation memory query (${VTYPE}$getri)',I0)")info
        call error(error_string)
      end if
    end if
    lwork = int(work2(1))

    allocate(work(lwork))
    call ${VTYPE}$getri(nn, aa, lda, ipiv, work, lwork, info)

    if (info < 0) then
      write (error_string, "('Failure in LU factorisation (${VTYPE}$getri), illegal argument at&
          & position ',I2)") info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
        write (error_string, "('Factor U is exactly zero in ${VTYPE}$getri, info flag is ',I0)")&
            & info
        call error(error_string)
      end if
    end if

  end subroutine getri_${NAME}$
#:endfor

#:for NAME, VTYPE, PREFIX, VPREC, CASE in [('real', 'real', 'sy', 'rsp', 's'),&
  & ('dreal', 'real', 'sy', 'rdp', 'd'), ('complex', 'complex', 'he', 'rsp', 'c'),&
  & ('dcomplex', 'complex', 'he', 'rdp', 'z')]
#:if PREFIX == 'sy'
  !> Computes the Bunch-Kaufman factorization of a symmetric matrix (${NAME}$).
#:else
  !> Computes the Bunch-Kaufman factorization of a hermitian matrix (${NAME}$).
#:endif
  subroutine ${PREFIX}$trf_${NAME}$(aa, ipiv, uplo, iError)

    !> Matrix to decompose, over-written on output
    ${VTYPE}$(${VPREC}$), intent(inout) :: aa(:,:)

    !> Interchanges of blocks on exit.
    integer, intent(out) :: ipiv(:)

    !> Signals whether upper (U) or lower (L) triangle should be used (default: lower).
    character, intent(in), optional :: uplo

    !> Info flag (0 = OK). If not set and an error occured, the subroutine stops.
    integer, intent(out), optional :: iError

    integer :: nn, info, lwork
    ${VTYPE}$(${VPREC}$), allocatable :: work(:)
    ${VTYPE}$(${VPREC}$) :: tmpwork(1)
    character :: uplo0

    if (present(uplo)) then
      uplo0 = uplo
    else
      uplo0 = "L"
    end if

    nn = size(aa, dim=2)
    lwork = -1
    call ${CASE}$${PREFIX}$trf(uplo0, nn, aa, nn, ipiv, tmpwork, lwork, info)
    if (info < 0) then
      write (error_string, "('Failure in factorisation memory query (${CASE}$${PREFIX}$trf),&
          & illegal argument at position ',I2)") info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
        write (error_string, "('Failure in factorisation memory query (${CASE}$${PREFIX}$trf)',&
            & I0)")info
        call error(error_string)
      end if
    end if
    lwork = int(tmpwork(1))
    allocate(work(lwork))
    call ${CASE}$${PREFIX}$trf(uplo0, nn, aa, nn, ipiv, work, lwork, info)
    if (info < 0) then
      write (error_string, "('Failure in factorisation memory query (${CASE}$${PREFIX}$trf),&
          & illegal argument at position ',I2)") info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
        write (error_string, "('Failure in factorisation memory query (${CASE}$${PREFIX}$trf)',&
            & I0)")info
        call error(error_string)
      end if
    end if

  end subroutine ${PREFIX}$trf_${NAME}$
#:endfor


  !> Computes the inverse of a symmetric matrix (real).
  subroutine sytri_real(aa, ipiv, uplo, info)

    !> Symmetric matrix to be inverted.
    real(rsp), intent(in) :: aa(:,:)

    !> Block interchanges as created by the sytrf() routine.
    integer, intent(in) :: ipiv(:)

    !> Upper ('U') or lower ('L') matrix (default: 'L')
    character, intent(in), optional :: uplo

    !> Info flag. If not present and an error occurs the subroutine stops.
    integer, intent(out), optional :: info

    integer :: info0, nn
    character :: uplo0
    real(rsp), allocatable :: work(:)

    if (present(uplo)) then
      uplo0 = uplo
    else
      uplo0 = "L"
    end if
    nn = size(aa, dim=1)
    allocate(work(max(1, 2 * nn)))
    call ssytri(uplo0, nn, aa, nn, ipiv, work, info0)
    if (present(info)) then
      info = info0
    elseif (info0 /= 0) then
      write(error_string, "(A,I10)") "Routine dsytri failed. Info: ", info0
      call error(error_string)
    end if

  end subroutine sytri_real


  !> Computes the inverse of a symmetric matrix (dreal).
  subroutine sytri_dreal(aa, ipiv, uplo, info)

    !> Symmetric matrix to be inverted.
    real(rdp), intent(in) :: aa(:,:)

    !> Block interchanges as created by the sytrf() routine.
    integer, intent(in) :: ipiv(:)

    !> upper or lower triangle
    character, intent(in), optional :: uplo

    !> Info flag. If not present and an error occurs the subroutine stops.
    integer, intent(out), optional :: info

    integer :: info0, nn
    character :: uplo0
    real(rdp), allocatable :: work(:)

    if (present(uplo)) then
      uplo0 = uplo
    else
      uplo0 = "L"
    end if
    nn = size(aa, dim=1)
    allocate(work(max(1, 2 * nn)))
    call dsytri(uplo0, nn, aa, nn, ipiv, work, info0)
    if (present(info)) then
      info = info0
    elseif (info0 /= 0) then
      write(error_string, "(A,I10)") "Routine dsytri failed. Info: ", info0
      call error(error_string)
    end if

  end subroutine sytri_dreal


  !> Computes the inverse of a Hermitian matrix (complex).
  subroutine hetri_complex(aa, ipiv, uplo, info)

    !> Symmetric matrix to be inverted.
    complex(rsp), intent(in) :: aa(:,:)

    !> Block interchanges as created by the sytrf() routine.
    integer, intent(in) :: ipiv(:)

    !> Upper ('U') or lower ('L') matrix (default: 'L')
    character, intent(in), optional :: uplo

    !> Info flag. If not present and an error occurs the subroutine stops.
    integer, intent(out), optional :: info

    integer :: info0, nn
    character :: uplo0
    complex(rsp), allocatable :: work(:)

    if (present(uplo)) then
      uplo0 = uplo
    else
      uplo0 = "L"
    end if
    nn = size(aa, dim=1)
    allocate(work(max(1, 2 * nn)))
    call chetri(uplo0, nn, aa, nn, ipiv, work, info0)
    if (present(info)) then
      info = info0
    elseif (info0 /= 0) then
      write(error_string, "(A,I10)") "Routine dsytri failed. Info: ", info0
      call error(error_string)
    end if

  end subroutine hetri_complex


  !> Computes the inverse of a Hermitian matrix (dreal).
  subroutine hetri_dcomplex(aa, ipiv, uplo, info)

    !> Hermitian matrix to be inverted.
    complex(rdp), intent(in) :: aa(:,:)

    !> Block interchanges as created by the sytrf() routine.
    integer, intent(in) :: ipiv(:)

    !> Upper ('U') or lower ('L') matrix (default: 'L')
    character, intent(in), optional :: uplo

    !> Info flag. If not present and an error occurs the subroutine stops.
    integer, intent(out), optional :: info

    integer :: info0, nn
    character :: uplo0
    complex(rdp), allocatable :: work(:)

    if (present(uplo)) then
      uplo0 = uplo
    else
      uplo0 = "L"
    end if
    nn = size(aa, dim=1)
    allocate(work(max(1, 2 * nn)))
    call zhetri(uplo0, nn, aa, nn, ipiv, work, info0)
    if (present(info)) then
      info = info0
    elseif (info0 /= 0) then
      write(error_string, "(A,I10)") "Routine dsytri failed. Info: ", info0
      call error(error_string)
    end if

  end subroutine hetri_dcomplex


  !> Single precision version of sytrs
  subroutine sytrs_real(A,B, nRow, uplo,iError)

    !> On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of
    !> A contains the upper triangular part of the matrix A, and the strictly lower triangular part
    !> of A is not referenced.  If UPLO = 'L', the leading N-by-N lower triangular part of A
    !> contains the lower triangular part of the matrix A, and the strictly upper triangular part of
    !> A is not referenced.  On exit, the block diagonal matrix D and the multipliers used to obtain
    !> the factor U or L
    real(rsp), intent(inout) :: A(:,:)

    !> On entry, the right hand side matrix B. On exit, the solution matrix X.
    real(rsp), intent(inout) :: B(:,:)

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of rows of
    !> the passed matrix)
    integer, intent(in), optional :: nRow

    !> upper or lower triangle of the matrix, defaults to lower
    character, intent(in), optional :: uplo

    !> Error flag. Zero on successfull exit. If not present, any lapack error causes program
    !> termination. If present, only fatal lapack errors with error flag < 0 cause abort.
    integer, intent(out), optional :: iError

    integer, allocatable :: ipiv(:)
    character :: iUplo
    integer :: nn, lda, ldb, info, lwork, nrhs
    real(rsp), allocatable :: work(:)
    real(rsp) :: work2(1)

    lda = size(A, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if

    ldb = size(b, dim=1)
    @:ASSERT(ldb >= nn)

    if (present(uplo)) then
      iUplo = uplo
    else
      iUplo = 'L'
    end if
    @:ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')

    @:ASSERT(size(A, dim=2) >= nn)
    nrhs = size(B, dim=2)

    allocate(ipiv(nn))

    lwork = -1
    call ssytrf(iUplo, nn, A, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))
    if (info == 0) then
      allocate(work(lwork))
      call ssytrf(iUplo, nn, A, lda, ipiv, work, lwork, info)
    end if

    if (info == 0) then
      call ssytrs(iUplo, nn, nrhs, A, lda, ipiv, B, ldb, info)
    end if

    if (present(iError)) then
      iError = info
    elseif (info /= 0) then
99130 format ('Solution failed because of error in sytrf or sytrs.',&
          & ' Info flag: ',i10)
      write (error_string, 99130) info
      call error(error_string)
    end if

  end subroutine sytrs_real


  !> Double precision version of sytrs
  subroutine sytrs_dble(A,B, nRow, uplo,iError)

    !> On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of
    !> A contains the upper triangular part of the matrix A, and the strictly lower triangular part
    !> of A is not referenced.  If UPLO = 'L', the leading N-by-N lower triangular part of A
    !> contains the lower triangular part of the matrix A, and the strictly upper triangular part of
    !> A is not referenced.  On exit, the block diagonal matrix D and the multipliers used to obtain
    !> the factor U or L
    real(rdp), intent(inout) :: A(:,:)

    !> On entry, the right hand side matrix B. On exit, the solution matrix X.
    real(rdp), intent(inout) :: B(:,:)

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of rows of
    !> the passed matrix)
    integer, intent(in), optional :: nRow

    !> upper or lower triangle of the matrix, defaults to lower
    character, intent(in), optional :: uplo

    !> Error flag. Zero on successfull exit. If not present, any lapack error causes program
    !> termination. If present, only fatal lapack errors with error flag < 0 cause abort.
    integer, intent(out), optional :: iError

    integer, allocatable :: ipiv(:)
    character :: iUplo
    integer :: nn, lda, ldb, info, lwork, nrhs
    real(rdp), allocatable :: work(:)
    real(rdp) :: work2(1)

    lda = size(A, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if

    ldb = size(b, dim=1)
    @:ASSERT(ldb >= nn)

    if (present(uplo)) then
      iUplo = uplo
    else
      iUplo = 'L'
    end if
    @:ASSERT(iUplo == 'u' .or. iUplo == 'U' .or. iUplo == 'l' .or. iUplo == 'L')

    @:ASSERT(size(A, dim=2) >= nn)
    nrhs = size(B, dim=2)

    allocate(ipiv(nn))

    lwork = -1
    call dsytrf(iUplo, nn, A, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))
    if (info == 0) then
      allocate(work(lwork))
      call dsytrf(iUplo, nn, A, lda, ipiv, work, lwork, info)
    end if

    if (info == 0) then
      call dsytrs(iUplo, nn, nrhs, A, lda, ipiv, B, ldb, info)
    end if

    if (present(iError)) then
      iError = info
    elseif (info /= 0) then
99130 format ('Solution failed because of error in sytrf or sytrs.',&
          & ' Info flag: ',i10)
      write (error_string, 99130) info
      call error(error_string)
    end if

  end subroutine sytrs_dble


  !> single precision version of larnv
  subroutine larnv_real(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> On entry, the seed of the random number generator; the array elements must be between 0 and
    !> 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    real(rsp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0
    call SLARNV( iDist, iSeed, n, x )
  end subroutine larnv_real


  !> double precision version of larnv
  subroutine larnv_dble(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> On entry, the seed of the random number generator; the array elements must be between 0 and
    !> 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    real(rdp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0d0
    call DLARNV( iDist, iSeed, n, x )
  end subroutine larnv_dble


  !> complex version of larnv
  subroutine larnv_cplx(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> On entry, the seed of the random number generator; the array elements must be between 0 and
    !> 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    complex(rsp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0
    call CLARNV( iDist, iSeed, n, x )
  end subroutine larnv_cplx


  !> double complex precision version of larnv
  subroutine larnv_dblecplx(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> INTEGER array, dimension (4) On entry, the seed of the random number generator; the array
    !> elements must be between 0 and 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    complex(rdp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0d0
    call ZLARNV( iDist, iSeed, n, x )
  end subroutine larnv_dblecplx


  !> real svd decomposition of matrix A into left and right vectors and singular values
  subroutine sgesvd_real(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    real(rsp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    real(rsp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rsp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    real(rsp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rsp), allocatable :: work(:)

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))

    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call sgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(work)

  end subroutine sgesvd_real


  !> double precision svd decomposition of matrix A into left and right vectors and singular values
  subroutine dgesvd_dble(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    real(rdp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    real(rdp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rdp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    real(rdp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rdp), allocatable :: work(:)

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))

    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call dgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(work)

  end subroutine dgesvd_dble

  !> complex svd decomposition of matrix A into left and right vectors and singular values
  subroutine cgesvd_cplx(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    complex(rsp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    complex(rsp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rsp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    complex(rsp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rsp), allocatable :: rwork(:)
    complex(rsp), allocatable :: work(:)

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = 2*min(m,n)+max(m,n)

    allocate(rwork(5*mn))
    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call cgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, rwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(rwork)
    deallocate(work)

  end subroutine cgesvd_cplx


  !> double complex svd decomposition of matrix A into left and right vectors and singular values
  subroutine zgesvd_dblecplx(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    complex(rdp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    complex(rdp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rdp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    complex(rdp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rdp), allocatable :: rwork(:)
    complex(rdp), allocatable :: work(:)

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = 2*min(m,n)+max(m,n)

    allocate(rwork(5*mn))
    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call zgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, rwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(rwork)
    deallocate(work)

  end subroutine zgesvd_dblecplx


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived routines requiring multiple LAPACK calls
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Inverts a matrix.
  subroutine matinv(aa, nRow, iError)

    !> Matrix to invert on entry, inverted matrix on exit
    real(dp), intent(inout) :: aa(:,:)

    !> Nr. of rows of the matrix (if different from size(aa, dim=1)
    integer, intent(in), optional :: nRow

    !> iError Error flag. Returns 0 on successfull operation. If this variable is not specified, any
    !> occuring error (e.g. singular matrix) stops the program.
    integer, intent(out), optional :: iError

    integer :: nn, info
    integer, allocatable :: ipiv(:)

    nn = size(aa, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= nn)
      nn = nRow
    end if
    @:ASSERT(size(aa, dim=2) >= nn)

    allocate(ipiv(nn))
    call getrf(aa, ipiv, nRow=nn, nColumn=nn, iError=iError)
    if (present(iError)) then
      if (iError /= 0) then
        write (error_string, "('Matrix inversion failed because of error in getrf. Info flag: ',&
            & I0)") info
        call warning(error_string)
      end if
    end if
    call getri(aa, ipiv, nRow=nn, iError=iError)
    if (present(iError)) then
      if (iError /= 0) then
        write (error_string, "('Matrix inversion failed because of error in getri. Info flag: ',&
            & I0)") info
        call warning(error_string)
      end if
    end if

  end subroutine matinv


#:for NAME, VTYPE, CASE, PREFIX in [('symmatinv', 'real', 'symmetric', 'sy'),&
  & ('hermatinv', 'complex', 'hermitian', 'he')]
  !> Inverts a ${CASE}$ matrix.
  subroutine ${NAME}$(aa, uplo, iError)

    !> Symmetric matrix to invert on entry, inverted matrix on exit.
    ${VTYPE}$(dp), intent(inout) :: aa(:,:)

    !> Upper ('U') or lower ('L') matrix. Default: 'L'.
    character, intent(in), optional :: uplo

    !> Info flag. If not specified and an error occurs, the subroutine will stop.
    integer, intent(out), optional :: iError

    integer :: nn
    integer, allocatable :: ipiv(:)

    nn = size(aa, dim=1)
    allocate(ipiv(nn))
    call ${PREFIX}$trf(aa, ipiv, uplo, iError=iError)
    if (present(iError)) then
      if (iError /= 0) then
        write (error_string, "('Matrix inversion failed because of error in ${PREFIX}$trf. Info&
            & flag: ', I0)") iError
        call warning(error_string)
      end if
    end if

    call ${PREFIX}$tri(aa, ipiv, uplo, info=iError)
    if (present(iError)) then
      if (iError /= 0) then
        write (error_string, "('Matrix inversion failed because of error in ${PREFIX}$tri. Info&
            & flag: ', I0)") iError
        call warning(error_string)
      end if
    end if

  end subroutine ${NAME}$
#:endfor

end module lapackroutines
