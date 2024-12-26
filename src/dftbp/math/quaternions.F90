!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines relating to quaternion rotation in 3D
module dftbp_math_quaternions
  use dftbp_common_accuracy, only : dp
  use dftbp_math_simplealgebra, only : cross3
  implicit none

  private
  public :: quaternionProduct, quaternionInvert, quaternionRotate, quaternionConstruct, rotate3
  public :: rotationMatrix

contains

  !> Product between two quaternions
  pure subroutine quaternionProduct(res,a,b)

    !> Result of product
    real(dp), intent(out) :: res(:)

    !> First quaternion
    real(dp), intent(in)  :: a(:)

    !> second quaternion
    real(dp), intent(in)  :: b(:)

    real(dp) axb(3)

    axb(:) = cross3(a(2:4), b(2:4))

    res(1) = a(1)*b(1) - dot_product(a(2:),b(2:))
    res(2:) = a(1)*b(2:) + b(1)*a(2:) + axb

  end subroutine quaternionProduct


  !> Form q^-1
  pure subroutine quaternionInvert(res)

    !> Initial quaternion, overwritten by result
    real(dp), intent(inout) :: res(:)

    real(dp) :: mag2

    mag2 = sum(res**2)

    res(2:) = -res(2:)
    res = res / mag2

  end subroutine quaternionInvert


  !> Rotates a 3 vector using a quaternion
  pure subroutine quaternionRotate(vec, quat)

    !> Initial 3 vector, overwritten by result
    real(dp), intent(inout) :: vec(:)

    !> Quaternion
    real(dp), intent(in)  :: quat(:)

    real(dp) :: qTmp(4,3)

    qTmp = 0.0_dp
    qTmp(2:,1) = vec

    qTmp(:,2) = quat
    call  quaternionInvert(qTmp(:,2))

    call quaternionProduct(qTmp(:,3),qTmp(:,1),qTmp(:,2))
    call quaternionProduct(qTmp(:,1),quat,qTmp(:,3))

    vec = qTmp(2:,1)

  end subroutine quaternionRotate


  !> Constructs a quaternion equivalent to rotation by an angle around an axis
  pure subroutine quaternionConstruct(res, angle, axis)

    !> The quaternion
    real(dp), intent(out) :: res(:)

    !> Angle of rotation
    real(dp), intent(in)  :: angle

    !> 3 axis of rotation, does not require normalisation
    real(dp), intent(in)  :: axis(:)

    real(dp) :: norm(3)

    norm = axis / sqrt(sum(axis**2))

    res(1) = cos(0.5_dp * angle)
    res(2:) = sin(0.5_dp * angle) * norm

  end subroutine quaternionConstruct


  !> Rotates a 3 vector by a specified angle around an axis
  pure subroutine rotate3(x, angle, axis)

    !> Vector, overwritten when rotated
    real(dp), intent(inout) :: x(:)

    !> Angle in radians
    real(dp), intent(in) :: angle

    !> Specified axis (does not require normalisation, but must be non-zero)
    real(dp), intent(in) :: axis(:)

    real(dp) :: quat(4)

    call quaternionConstruct(quat,angle,axis)
    call quaternionRotate(x,quat)

  end subroutine rotate3


  !> 3D rotation matrix from a quaternion (Euler-Rodrigues formula)
  pure function rotationMatrix(q) result(r)

    !> Quaternion
    real(dp), intent(in) :: q(0:3)

    !> Resulting rotation matrix
    real(dp) :: r(3,3)

    integer :: ii

    do ii = 1, 3
      r(ii,ii) = q(0)**2 + q(ii)**2
    end do
    r(2,1) = q(1)*q(2)-q(0)*q(3)
    r(3,1) = q(1)*q(3)+q(0)*q(2)
    r(3,2) = q(2)*q(3)-q(0)*q(1)
    r(1,2) = q(1)*q(2)+q(0)*q(3)
    r(1,3) = q(1)*q(3)-q(0)*q(2)
    r(2,3) = q(2)*q(3)+q(0)*q(1)
    r(:,:) = 2.0_dp * r
    do ii = 1, 3
      r(ii,ii) = r(ii,ii) - 1.0_dp
    end do

  end function rotationMatrix

end module dftbp_math_quaternions
