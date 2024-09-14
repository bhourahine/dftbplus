!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Higer accuracy summation of values
!! TODO: add pairwise summation algorithm (OMP and also use Kahan for base sum at end)
module dftbp_math_summation
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: KahanBabushkaKleinSum, dot

contains

  !> Dot product
  pure function dot(a,b)

    !> Vector to product
    real(dp), intent(in) :: a(:)

    !> Vector to product
    real(dp), intent(in) :: b(:)

    !> Resulting product
    real(dp) :: dot

    real(dp) :: c(size(a))

    c(:) = a * b
    dot = KahanBabushkaKleinSum(c)

  end function dot


  !> Higher order Kahan summation, should have O(eps) error. Note, -Ofast dangerous, but it seems OK
  !! up to -O3 with gcc 13.2.
  !! Klein "A generalized Kahan–Babuška-Summation-Algorithm" Computing 76, 279–293
  !! (2006). doi:10.1007/s00607-005-0139-x
  pure function KahanBabushkaKleinSum(array)

    !> Array of values to sum
    real(dp), intent(in) :: array(:)

    !> Resulting sum
    real(dp) :: KahanBabushkaKleinSum

    real(dp) :: sum, cs, ccs, c, cc, t
    integer :: ii

    sum = 0.0
    cs  = 0.0
    ccs = 0.0
    c   = 0.0
    cc  = 0.0

    do ii = 1, size(array)
      t = sum + array(ii)
      if (abs(sum) >= abs(array(ii))) then
        c = (sum - t) + array(ii)
      else
        c = (array(ii) - t) + sum
      endif
      sum = t
      t = cs + c
      if (abs(cs) >= abs(c)) then
        cc = (cs - t) + c
      else
        cc = (c - t) + cs
      endif
      cs = t
      ccs = ccs + cc
    end do

    KahanBabushkaKleinSum = sum + cs + ccs

  end function KahanBabushkaKleinSum

end module dftbp_math_summation
