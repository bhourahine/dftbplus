!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains routines to locate a value in a sorted array using binary search
module dftbp_math_binarysearch
  use dftbp_common_accuracy, only : dp
  use dftbp_common_status, only : TStatus
  use dftbp_math_comparisions, only : eq, gt, ge, lt, le
  use dftbp_io_message, only : error
  implicit none

  private
  public :: search_int, search_asc_real_geq, search_asc_real_gt, search_des_real_geq,&
      & search_des_real_gt, isPresent, first, last, firstGreater, lastLesser


  interface isPresent
    module procedure :: isPresent_intScalar
    module procedure :: isPresent_intVector
  end interface isPresent

  interface first
    module procedure :: first_intScalar
    module procedure :: first_intVector
  end interface first


  interface last
    module procedure :: last_intScalar
    module procedure :: last_intVector
  end interface last


  interface firstGreater
    module procedure :: firstGreater_intScalar
    module procedure :: firstGreater_intVector
  end interface firstGreater


  interface lastLesser
    module procedure :: lastLesser_intScalar
    module procedure :: lastLesser_intVector
  end interface lastLesser


contains


  !> Integer case for binary search of sorted values to find the jj such that
  !! xVal < xx(jj+1), i.e., the last occurance of the value in an array, or if not present in the
  !! array, the last element smaller than xVal. If xVal < xx(1), jj = 0
  subroutine search_int(jj, xx, xVal)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    integer, intent(in) :: xx(:)

    !> Value to locate jj for
    integer, intent(in) :: xVal

    integer :: jlower, jupper, jcurr

    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if (xx(jcurr) <= xVal) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_int


  !======================================
  ! Ascending ordered real array versions
  !======================================


  !> Real case for binary search of ascending sorted values to find the jj such that
  !! xx(jj) + tol <= xVal, i.e., the last element equal to or smaller than xVal. If xVal < xx(1) -
  !! tol, jj = 0
  subroutine search_asc_real_geq(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if (xVal - xx(jcurr) >= -tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_asc_real_geq


  !> Real case for binary search of ascending sorted values to find the jj such that
  !! xVal < xx(jj+1) - tol, i.e. the last occurance of a value smaller than xVal. If xVal < xx(1) -
  !! tol, jj = 0
  subroutine search_asc_real_gt(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if (xVal - xx(jcurr) > tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_asc_real_gt


  !=======================================
  ! Descending ordered real array versions
  !=======================================


  !> Real case for binary search of decending sorted values to find the first jj such that
  !! xx(jj) + tol >= xVal
  subroutine search_des_real_geq(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if ((xx(jcurr) - xVal) >= -tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_des_real_geq


  !> Real case for binary search of descending sorted values to find first element jj such that
  !! xx(jj) - tol > xVal
  subroutine search_des_real_gt(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if ((xx(jcurr) - xVal) > tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_des_real_gt


  !> Tests if a value is present in an ascending sorted array
  pure function isPresent_intScalar(low, high, xx, xVal) result(isPresent)
    integer, intent(in) :: low, high, xx(:), xVal
    logical :: isPresent
    integer :: jCurr, jCurrVal, jLower, jUpper
    jLower = low
    jUpper = high
    isPresent = .false.
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower) / 2
      jCurrVal = xx(jCurr)
      if (jCurrVal < xVal) then
        jLower = jCurr + 1
      elseif (jCurrVal > xVal) then
        jUpper = jCurr - 1
      else
        isPresent = .true.
        return
      end if
    end do
  end function isPresent_intScalar


  !> Tests if a value is present in an ascending sorted array
  pure function isPresent_intVector(low, high, xx, xVal, fields) result(isPresent)
    integer, intent(in) :: low, high, xx(:,:), xVal(:), fields(:)
    logical :: isPresent
    integer :: jCurr, jLower, jUpper
    integer, allocatable :: jCurrVal(:)
    jLower = low
    jUpper = high
    isPresent = .false.
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower) / 2
      jCurrVal = xx(:, jCurr)
      if (lt(jCurrVal, xVal, fields)) then
        jLower = jCurr + 1
      else if (gt(jCurrVal, xVal, fields)) then
        jUpper = jCurr - 1
      else
        isPresent = .true.
        return
      end if
    end do
  end function isPresent_intVector


  !> Find first occurrence index of xVal in an ascending sorted array
  pure function first_intScalar(low, high, xx, xVal) result(jj)
    integer, intent(in) :: low, high, xx(:), xVal
    integer :: jj
    integer :: jCurr, jCurrVal, jLower, jUpper
    jLower = low
    jUpper = high
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(jCurr)
      if (jCurrVal < xVal) then
        jLower = jCurr + 1
      elseif (jCurrVal > xVal) then
        jUpper = jCurr - 1
      else
        jj = jCurr
        ! look lower than current best match
        jUpper = jCurr - 1
      end if
    end do
  end function first_intScalar


  !> Find first occurrence index of xVal in an ascending sorted array
  pure function first_intVector(low, high, xx, xVal, fields) result(jj)
    integer, intent(in) :: low, high, xx(:,:), xVal(:), fields(:)
    integer :: jj
    integer :: jCurr, jLower, jUpper
    integer, allocatable :: jCurrVal(:)
    jLower = low
    jUpper = high
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(:, jCurr)
      if (lt(jCurrVal, xVal, fields)) then
        jLower = jCurr + 1
      elseif (gt(jCurrVal, xVal, fields)) then
        jUpper = jCurr - 1
      else
        jj = jCurr
        ! look lower than current best match
        jUpper = jCurr - 1
      end if
    end do
  end function first_intVector


  !> Find last occurrence index of xVal in xx
  pure function last_intScalar(low, high, xx, xVal) result(jj)
    integer, intent(in) :: low, high, xx(:), xVal
    integer :: jj
    integer :: jCurr, jCurrVal, jLower, jUpper
    jLower = low
    jUpper = high
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(jCurr)
      if (jCurrVal < xVal) then
        ! look higher
        jLower = jCurr + 1
      elseif (jCurrVal > xVal) then
        ! look lower
        jUpper = jCurr - 1
      else
        jj = jCurr
        ! look higher than current best match
        jLower = jCurr + 1
      end if
    end do
  end function last_intScalar


  !> Find last occurrence index of xVal in xx
  pure function last_intVector(low, high, xx, xVal, fields) result(jj)
    integer, intent(in) :: low, high, xx(:,:), xVal(:), fields(:)
    integer :: jj
    integer :: jCurr, jLower, jUpper
    integer, allocatable :: jCurrVal(:)
    jLower = low
    jUpper = high
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(:,jCurr)
      if (lt(jCurrVal, xVal, fields)) then
        ! look higher
        jLower = jCurr + 1
      elseif (gt(jCurrVal, xVal, fields)) then
        ! look lower
        jUpper = jCurr - 1
      else
        jj = jCurr
        ! look higher than current best match
        jLower = jCurr + 1
      end if
    end do
  end function last_intVector


  !> Find index of first occurrence of element greater than xVal in xx
  pure function firstGreater_intScalar(low, high, xx, xVal) result(jj)
    integer, intent(in) :: low, high, xx(:), xVal
    integer :: jj
    integer :: jCurr, jCurrVal, jLower, jUpper
    jLower = low
    jUpper = high
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(jCurr)
      if (jCurrVal <= xVal) then
        ! look higher
        jLower = jCurr + 1
      elseif (jCurrVal > xVal)  then
        jj = jCurr
        ! look lower than current best match
        jUpper = jCurr - 1
      end if
    end do
  end function firstGreater_intScalar


  !> Find index of first occurrence of element greater than xVal in xx
  pure function firstGreater_intVector(low, high, xx, xVal, fields) result(jj)
    integer, intent(in) :: low, high, xx(:,:), xVal(:), fields(:)
    integer :: jj
    integer :: jCurr, jLower, jUpper
    integer, allocatable :: jCurrVal(:)
    jLower = low
    jUpper = high
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(:, jCurr)
      if (le(jCurrVal, xVal, fields)) then
        ! look higher
        jLower = jCurr + 1
      elseif (gt(jCurrVal, xVal, fields))  then
        jj = jCurr
        ! look lower than current best match
        jUpper = jCurr - 1
      end if
    end do
  end function firstGreater_intVector


  !> Find index of last occurrence of element less than xVal in xx
  pure function lastLesser_intScalar(low, high, xx, xVal) result(jj)
    integer, intent(in) :: low, high, xx(:), xVal
    integer :: jj
    integer :: jCurr, jCurrVal, jLower, jUpper
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(jCurr)
      if (jCurrVal < xVal) then
        jj = jCurr
        ! look higher than current best match
        jLower = jCurr + 1
      elseif (jCurrVal >= xVal) then
        ! look lower
        jUpper = jCurr - 1
      end if
    end do
  end function lastLesser_intScalar


  !> Find index of last occurrence of element less than xVal in xx
  pure function lastLesser_intVector(low, high, xx, xVal, fields) result(jj)
    integer, intent(in) :: low, high, xx(:,:), xVal(:), fields(:)
    integer :: jj
    integer :: jCurr, jLower, jUpper
    integer, allocatable :: jCurrVal(:)
    jj = -1
    do while (jLower <= jUpper)
      jCurr = jLower + (jUpper - jLower + 1) / 2
      jCurrVal = xx(:, jCurr)
      if (lt(jCurrVal, xVal, fields)) then
        jj = jCurr
        ! look higher than current best match
        jLower = jCurr + 1
      elseif (ge(jCurrVal, xVal, fields)) then
        ! look lower
        jUpper = jCurr - 1
      end if
    end do
  end function lastLesser_intVector


end module dftbp_math_binarysearch
