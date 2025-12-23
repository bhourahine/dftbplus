!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_math_angmomentum
  use dftbp_common_accuracy, only : dp
  use dftbp_math_angmomentum, only : coupledBasis, clebschGordan
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("cgCoefficients")

    integer, parameter :: jmax = 3
    integer :: j, m
    real(dp) :: cg, analytic

    write(*,*)
    write(*,*)'Symmetries'
    do j = 0, jmax
      do m = -j, j
        cg = clebschGordan(real(j,dp),real(m,dp),real(j,dp),real(-m,dp),0.0_dp,0.0_dp)
        analytic = (-1.0_dp**(j-m)) / sqrt(2.0_dp*real(j,dp)+1.0_dp)
        write(*,*)j,m,cg, analytic
        @:ASSERT(abs(cg - analytic) < epsilon(0.0_dp))
      end do
    end do
    #!@:ASSERT(.false.)

  $:END_TEST()


  $:TEST("coupledBasis")

    integer, parameter :: lmax = 3
    real(dp) :: matReal(2*lmax+1, 2*lmax+1)
    integer :: l

    do l = 0, 3
      write(*,"(1X,A,I0)")'L=',l
      call coupledBasis(matReal, l)
    end do

    #!@:ASSERT(.false.)

  $:END_TEST()



  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("angmomentum", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_angmomentum
