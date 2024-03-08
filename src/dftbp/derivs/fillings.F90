!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for determining filled and empty bands
module dftbp_derivs_fillings
  use dftbp_common_accuracy, only : dp
  implicit none

  public

contains

  pure function nIndependentHam(nSpin)

    integer, intent(in) :: nSpin
    integer :: nIndependentHam

    select case(nSpin)
    case(1,4)
      nIndependentHam = 1
    case(2)
      nIndependentHam = 2
    end select

  end function nIndependentHam


  pure function maximumFillings(nSpin)

    integer, intent(in) :: nSpin
    real(dp) :: maximumFillings

    select case(nSpin)
    case(1)
      maximumFillings = 2.0_dp
    case(2,4)
      maximumFillings = 1.0_dp
    end select

  end function maximumFillings


  subroutine filledOrEmpty(nFilled, nEmpty, nIndepHam, nKpts, filling, nOrbs, maxFill)

    integer, allocatable, intent(out) :: nFilled(:,:), nEmpty(:,:)

    integer, intent(in) :: nIndepHam, nKpts, nOrbs
    real(dp), intent(in) :: filling(:,:,:), maxFill

    integer :: ik, iLev, iS

    allocate(nFilled(nIndepHam, nKpts))
    allocate(nEmpty(nIndepHam, nKpts))

    nFilled(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKPts
        do iLev = 1, nOrbs
          if ( filling(iLev, iK, iS) < epsilon(1.0) ) then
            ! assumes Fermi filling, so above this is empty
            nFilled(iS, iK) = iLev - 1
            exit
          end if
        end do
        ! check if channel is fully filled
        if (nFilled(iS, iK) < 0) then
          nFilled(iS, iK) = nOrbs
        end if
      end do
    end do
    nEmpty(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKpts
        do iLev = 1, nOrbs
          if ( abs( filling(iLev, iK, iS) - maxFill ) > epsilon(1.0)) then
            ! assumes Fermi filling, so this is filled
            nEmpty(iS, iK) = iLev
            exit
          end if
        end do
        !> Check is channel is empty
        if (nEmpty(iS, iK) < 0) then
          nEmpty(iS, iK) = 1
        end if
      end do
    end do

  end subroutine filledOrEmpty

end module dftbp_derivs_fillings
