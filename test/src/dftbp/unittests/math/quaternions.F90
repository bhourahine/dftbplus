!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("quaternions")
  use dftbp_common_accuracy, only : dp
  use dftbp_math_quaternions, only : rotationMatrix, quaternionConstruct
  implicit none

#:contains

#:block TEST_FIXTURE("rotate")

     real(dp) :: q(4), r(3,3), v(3)

  #:contains

    #:block TEST("identity")
       v(:) = [0.5_dp, -.75_dp, 1.0_dp]
       call quaternionConstruct(q, 0.123456_dp, v)
       r(:,:) = rotationMatrix(q)
       @:ASSERT(all(abs(v - matmul(r, v)) < 128_dp * epsilon(0.0_dp)))
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
