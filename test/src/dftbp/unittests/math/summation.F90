!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("summation")
  use dftbp_common_accuracy, only : dp, rdp, rsp
  use dftbp_math_summation, only : KahanBabushkaKleinSum
  implicit none

#:contains

  #:block TEST_FIXTURE("underflow")

  real(dp) :: vector(6) = [1.0_dp,0.499_dp*huge(1.0_rdp),1.0_dp*huge(1.0_rsp),1.0_dp,&
      &-0.499_dp*huge(1.0_rdp),-1.0_dp*huge(1.0_rsp)]

  #:contains

    #:block TEST("sum")
       @:ASSERT(.true.)
       #!@:ASSERT(abs(KahanBabushkaKleinSum(vector) - 2.0_dp) < epsilon(1.0_dp))
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
