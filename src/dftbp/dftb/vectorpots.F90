!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module for vector potentials
module dftbp_dftb_vectorpots
  use dftbp_common_accuracy, only : dp
  use dftbp_math_simplealgebra, only : cross3
  implicit none

  private
  public :: circular, xyPlane

contains

  !> Circular (symmetric) Couomb gauge for arbitrary oriented B field and origin
  pure function circular(B, r, O) result(A)

    !> z aligned magnetic field
    real(dp), intent(in) :: B(3)

    !> Position
    real(dp), intent(in) :: r(3)

    !> Gauge origin
    real(dp), intent(in) :: O(3)

    !> Vector potential
    real(dp) :: A(3)

    A = 0.5_dp * cross3(B, r - O)

  end function circular


  !> Fairly general range of Coulomb gauges for the vector potential in the xy plane for Bz field,
  !> dependinging on choice of beta
  pure function xyPlane(Bz, x, beta) result(A)

    !> z aligned magnetic field
    real(dp), intent(in) :: Bz

    !> Position
    real(dp), intent(in) :: x(3)

    !> beta = 0 circular (symmetric) gauge, = +/-1 Landau gauge with A aligned x and y respectively
    real(dp), intent(in) :: beta

    !> Vector potential
    real(dp) :: A(3)

    A(:) = Bz * 0.5_dp * [-x(2) * (1.0_dp + beta), x(1) * (1.0_dp - beta), 0.0_dp]

  end function xyPlane

end module dftbp_dftb_vectorpots
