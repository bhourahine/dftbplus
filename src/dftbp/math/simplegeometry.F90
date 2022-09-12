!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Simple geometry calculations
module dftbp_math_simplegeometry
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: dist2Segment

contains

  !> Distance of a point to a line segment between points p1 and p2
  pure function dist2Segment(p1, p2, q) result(dist)

    !> Start of segment
    real(dp), intent(in) :: p1(:)

    !> End of segment
    real(dp), intent(in) :: p2(:)

    !> Point to measure distance from
    real(dp), intent(in) :: q(:)

    !> Resulting distance to segment
    real(dp) :: dist

    real(dp) :: u(size(p1)), v(size(p1)), p(size(p1)), d

    u(:) = p2 - p1
    v(:) = q - p1

    d = dot_product(u, v) / dot_product(u,u) ! projection onto line segment
    if (d < 0.0_dp) then ! projection is before start of segment
      dist = sqrt(dot_product(q-p1,q-p1))
    else if (d > 1.0_dp) then ! projection is after end of seqment
      dist = sqrt(dot_product(q-p2,q-p2))
    else ! projection is inside segment
      p(:) = p1 + d * u ! nearest point inside segment
      dist = sqrt(dot_product(q-p,q-p))
    end if

  end function dist2Segment

end module dftbp_math_simplegeometry
