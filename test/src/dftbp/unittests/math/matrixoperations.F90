!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("matrixoperations")
  use dftbp_common_accuracy, only : dp
  use dftbp_math_matrixoperations, only : triangleCopySquareMatrix, blockTriangleCopyHS
  implicit none

#:contains

  #:block TEST_FIXTURE("denseSymmetry")

     integer :: ii, jj, kk
     integer, parameter :: n = 3
     real(dp) :: matReal(n,n)
     complex(dp) :: matCplx(n,n)

  #:contains

    #:block TEST("real")
       @:ASSERT(size(matReal, dim=1) == size(matReal, dim=2))
       matReal(:,:) = 0.0_dp
       kk = 0
       do ii = 1, size(matReal, dim=1)
         do jj = ii, size(matReal, dim=2)
           kk = kk + 1
           matReal(jj, ii) = real(kk,dp)
         end do
       end do
       call triangleCopySquareMatrix(matReal)
       @:ASSERT(all(matReal > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matReal - transpose(matReal)) < epsilon(0.0_dp)))
    #:endblock

    #:block TEST("complex")
       @:ASSERT(size(matCplx, dim=1) == size(matCplx, dim=2))
       matCplx(:,:) = cmplx(0,0,dp)
       kk = 0
       do ii = 1, size(matCplx, dim=1)
         do jj = ii, size(matCplx, dim=2)
           kk = kk + 1
           matCplx(jj, ii) = real(kk,dp)
         end do
       end do
       kk = 0
       do ii = 1, size(matCplx, dim=1)
         do jj = ii+1, size(matCplx, dim=2)
           kk = kk + 1
           matCplx(jj, ii) = matCplx(jj, ii) + cmplx(0,kk,dp)
         end do
       end do
       call triangleCopySquareMatrix(matCplx)
       @:ASSERT(all(abs(matCplx) > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matCplx - transpose(conjg(matCplx))) < epsilon(0.0_dp)))
    #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("blockSymmetry")

     integer :: ii, jj, kk
     integer, parameter :: n = 5
     integer, parameter :: nAt = 2
     real(dp) :: matReal(n,n)
     complex(dp) :: matCplx(n,n), matPauli(2*n,2*n)
     integer, parameter :: iAtomStart(nAt+1) = [1, 2, n+1]

  #:contains

    #:block TEST("real")
       @:ASSERT(size(matReal, dim=1) == size(matReal, dim=2))
       matReal(:,:) = 0.0_dp
       kk = 0
       do ii = 1, size(matReal, dim=1)
         do jj = ii, size(matReal, dim=2)
           kk = kk + 1
           matReal(jj, ii) = real(kk,dp)
         end do
       end do
       ! symmetrize diagonal atomic blocks, as call being tested does not touch these
       do ii = 1, nAt
         do jj = iAtomStart(ii), iAtomStart(ii+1)-1
           do kk = jj+1, iAtomStart(ii+1)-1
             matReal(jj,kk) = matReal(kk,jj)
           end do
         end do
       end do
       call blockTriangleCopyHS(matReal, iAtomStart)
       @:ASSERT(all(matReal > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matReal - transpose(matReal)) < epsilon(0.0_dp)))
    #:endblock


    #:block TEST("complex")
       @:ASSERT(size(matCplx, dim=1) == size(matCplx, dim=2))
       matCplx(:,:) = cmplx(0,0,dp)
       kk = 0
       do ii = 1, size(matCplx, dim=1)
         do jj = ii, size(matCplx, dim=2)
           kk = kk + 1
           matCplx(jj, ii) = real(kk,dp)
         end do
       end do
       kk = 0
       do ii = 1, size(matCplx, dim=1)
         do jj = ii+1, size(matCplx, dim=2)
           kk = kk + 1
           matCplx(jj, ii) = matCplx(jj, ii) + cmplx(0,kk,dp)
         end do
       end do
       ! Hermitian symmetrize diagonal atomic blocks, as call being tested does not touch these
       do ii = 1, nAt
         do jj = iAtomStart(ii), iAtomStart(ii+1)-1
           do kk = jj+1, iAtomStart(ii+1)-1
             matCplx(jj,kk) = conjg(matCplx(kk,jj))
           end do
         end do
       end do
       call blockTriangleCopyHS(matCplx, iAtomStart)
       @:ASSERT(all(abs(matCplx) > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matCplx - transpose(conjg(matCplx))) < epsilon(0.0_dp)))
    #:endblock


    #:block TEST("pauli")
       @:ASSERT(size(matPauli, dim=1) == size(matPauli, dim=2))
       matPauli(:,:) = cmplx(0,0,dp)
       kk = 0
       do ii = 1, size(matPauli, dim=1)
         do jj = ii, size(matPauli, dim=2)
           kk = kk + 1
           matPauli(jj, ii) = real(kk,dp)
         end do
       end do
       kk = 0
       do ii = 1, size(matPauli, dim=1)
         do jj = ii+1, size(matPauli, dim=2)
           kk = kk + 1
           matPauli(jj, ii) = matPauli(jj, ii) + cmplx(0,kk,dp)
         end do
       end do
       ! Hermitian symmetrize diagonal atomic blocks, as call being tested does not touch these
       do ii = 1, nAt
         do jj = iAtomStart(ii), iAtomStart(ii+1)-1
           do kk = jj+1, iAtomStart(ii+1)-1
             matPauli(jj,kk) = conjg(matPauli(kk,jj))
             matPauli(n+jj,n+kk) = conjg(matPauli(n+kk,n+jj))
           end do
         end do
       end do
       call blockTriangleCopyHS(matPauli, iAtomStart)
       @:ASSERT(all(abs(matPauli) > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matPauli - transpose(conjg(matPauli))) < epsilon(0.0_dp)))
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
