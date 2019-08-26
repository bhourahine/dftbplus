!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> contains some miscellaneous quantum mechanics related bits and pieces.
module dftbp_qm
  use dftbp_assert
  use dftbp_accuracy, only : dp, rsp
  use dftbp_blasroutines, only : gemm
  use dftbp_message, only : error
#:if WITH_SCALAPACK
  use dftbp_mpifx
  use dftbp_scalapackfx
#:endif
  implicit none

  private
  public :: makeSimiliarityTrans, commutator, nonOrthogCommutator
#!:if WITH_SCALAPACK
#!  public :: makeSimiliarityTrans_pblas
#!:endif

#:set FLAVOURS = [('Real', 'real'), ('Cmplx', 'complex')]

  !> perform a similarity (or unitary) transformation of a matrix X := U X U^dag, where U and U^dag
  !> can be interchanged if neccessary
  interface makeSimiliarityTrans
  #:for NAME, _ in FLAVOURS
    module procedure U_${NAME}$
  #:endfor
  end interface makeSimiliarityTrans

  !> Evaluate a commutator between matrices [A,B]
  interface commutator
  #:for NAME, _ in FLAVOURS
    module procedure commute${NAME}$
  #:endfor
  end interface commutator

  !> Evaluate a non-orthogonal commutator between matrices [FP,S]
  interface nonOrthogCommutator
  #:for NAME, _ in FLAVOURS
    module procedure commuteFPS${NAME}$
  #:endfor
  end interface nonOrthogCommutator


contains

#:for NAME, VAR in FLAVOURS

  !> Unitary transformation of a matrix matrix X := U X U^dag or X := U^dag X U for sides L or R
  !> respectively
  subroutine U_${NAME}$(xx, uu, side)

    !> matrix in original basis, overwritten on return.
    ${VAR}$(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    ${VAR}$(dp), intent(in) :: uu(:,:)

    !> which transform order to use, i.e. to which side of the matrix the original unitary is
    !> applied
    character(1), intent(in), optional :: side

    ${VAR}$(dp) :: work(size(xx,dim=1), size(xx,dim=2))
    character(1) :: iSide

  #:if VAR == 'real'
    character(1), parameter :: op = 'T'
  #:else
    character(1), parameter :: op = 'C'
  #:endif

    if (present(side)) then
      iSide(:) = side
    else
      iSide(:) = 'L'
    end if

    select case(iSide)

    case ('L', 'l')

      call gemm(work, xx, uu, transB = op)
      call gemm(xx, uu, work)

    case('R', 'r')

      call gemm(work, xx, uu)
      call gemm(xx, uu, work, transA = op)

    case default

      call error("Unknown side in U_${NAME}$ : "//trim(iSide))

    end select

  end subroutine U_${NAME}$


  !> Evaluate the commutator [A,B] for ${VAR}$ matrices
  subroutine commute${NAME}$(C, A, B)

    !> Resulting matrix
    ${VAR}$(dp), intent(out) :: C(:,:)

    !> First matrix
    ${VAR}$(dp), intent(in) :: A(:,:)

    !> Second matrix
    ${VAR}$(dp), intent(in) :: B(:,:)

  #:if VAR == 'real'
    !> value of +1
    real(dp), parameter :: one = 1.0_dp

    !> value of -1
    real(dp), parameter :: minusOne = -1.0_dp
  #:else
    !> value of +1
    complex(dp), parameter :: one = (1.0_dp,0.0_dp)

    !> value of -1
    complex(dp), parameter :: minusOne = (-1.0_dp, 0.0_dp)
  #:endif

    !C = matmul(A,B) - matmul(B,A)
    call gemm(C, A, B)
    call gemm(C, B, A, alpha = minusOne, beta = one)

  end subroutine commute${NAME}$


  !> Evaluate the Pulay commutator FPS-SPF for ${VAR}$ matrices, being a direct error measure for
  !> self-consistency with idempotent density matrices, or if S^-1 is passed the time derivative
  subroutine commuteFPS${NAME}$(e, F, P, S)

    !> Resulting matrix
    ${VAR}$(dp), intent(out) :: e(:,:)

    !> Hamiltonian (Fock) matrix
    ${VAR}$(dp), intent(in) :: F(:,:)

    !> Density matrix
    ${VAR}$(dp), intent(in) :: P(:,:)

    !> Overlap matrix
    ${VAR}$(dp), intent(in) :: S(:,:)

    ${VAR}$(dp), allocatable :: work(:,:)

  #:if VAR == 'real'
    character(1), parameter :: op = 'T'
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: minusOne = -1.0_dp
  #:else
    character(1), parameter :: op = 'C'
    complex(dp), parameter :: one = (1.0_dp,0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp,0.0_dp)
    complex(dp), parameter :: minusOne = (-1.0_dp, 0.0_dp)
  #:endif

    allocate(work(size(s,dim=1), size(s,dim=1)))

    e(:,:) = zero
    call gemm(work, S, F)
    call gemm(e, work, P)
    call gemm(e, P, work, alpha = minusOne, beta = one, transB = op)

    ! e(:,:) = matmul(F,matmul(P,S)) - matmul(S,matmul(P,F))

  end subroutine commuteFPS${NAME}$

#:endfor

end module dftbp_qm
