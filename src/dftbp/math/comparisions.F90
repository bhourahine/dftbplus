!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Various types of comparison routines
module dftbp_math_comparisions
  implicit none

  private
  public :: eq, gt, ge, lt, le

contains

  !> Integer multiple-key .eq.
  pure function eq(x, y, keys)

    !> First values
    integer, intent(in) :: x(:)

    !> Second values
    integer, intent(in) :: y(:)

    !> Which keys to use
    integer, intent(in) :: keys(:)

    logical :: eq

    eq = all(x(keys) == y(keys))

  end function eq


  !> Integer multiple-key .gt.
  pure function gt(x, y, keys)

    integer, intent(in) :: x(:)
    integer, intent(in) :: y(:)
    integer, intent(in) :: keys(:)

    logical :: gt

    integer :: ii, jj

    gt = .false.
    do ii = 1, size(keys)
      jj = keys(ii)
      if (x(jj) > y(jj)) then
        gt = .true.
        return
      else if (x(jj) < y(jj)) then
        gt = .false.
        return
      end if
    end do

  end function gt


  !> Integer multiple-key .ge.
  pure function ge(x, y, keys)

    integer, intent(in) :: x(:)
    integer, intent(in) :: y(:)
    integer, intent(in) :: keys(:)

    logical :: ge

    ge = gt(x, y, keys) .or. eq(x, y, keys)

  end function ge


  !> Integer multiple-key .lt.
  pure function lt(x, y, keys)

    integer, intent(in) :: x(:)
    integer, intent(in) :: y(:)
    integer, intent(in) :: keys(:)

    logical :: lt

    integer :: ii, jj

    lt = .false.
    do ii = 1, size(keys)
      jj = keys(ii)
      if (x(jj) < y(jj)) then
        lt = .true.
        return
      else if (x(jj) > y(jj)) then
        lt = .false.
        return
      end if
    end do

  end function lt


  !> Integer multiple-key .le.
  pure function le(x, y, keys)

    integer, intent(in) :: x(:)
    integer, intent(in) :: y(:)
    integer, intent(in) :: keys(:)

    logical :: le

    le = lt(x, y, keys) .or. eq(x, y, keys)

  end function le


end module dftbp_math_comparisions
