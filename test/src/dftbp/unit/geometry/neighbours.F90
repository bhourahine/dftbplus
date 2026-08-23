!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"
module test_geometry_neighbours
  use dftbp_common_accuracy, only : dp
  use dftbp_geometry_neighbours, only : TNeighborFinder
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  $:TEST("speedtest")
    integer, parameter :: nAt = 10000
    real(dp), allocatable :: coords(:,:)

    @:ASSERT(.true.)

  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("geometry", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_geometry_neighbours
