!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines relating to wavefunction gauge
module dftbp_math_gauge
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag, pi
  use dftbp_type_densedescr, only : TDenseDescr
  implicit none

  private
  public :: coordinateIncluding

  !> Transform eigenvectors to a better convention to evaluate properties like optical matrix
  !! elements or Berry phases
  interface coordinateIncluding
    module procedure coordinateIncluding_serial
  end interface coordinateIncluding

contains

  !> Serial transformation routines
  subroutine coordinateIncluding_serial(eigVecs, kPoint, nAtoms, denseDesc, coords)

    !> Eigenvectors in cell translation gauge
    complex(dp), intent(inout) :: eigVecs(:,:)

    !> Relative coordinates of the current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Number of central cell atoms
    integer, intent(in) :: nAtoms

    !> Descriptors for dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> Coordinates of the central cell atoms
    real(dp), intent(in) :: coords(:,:)

    complex(dp) :: phase
    integer :: iAt, iEnd, iStart
    real(dp) :: kVec(3)

    kVec(:) = 2.0_dp * pi * kPoint
    do iAt = 1, nAtoms
      iStart = denseDesc%iAtomStart(iAt)
      iEnd = denseDesc%iAtomStart(iAt+1)-1
      phase = exp(imag * dot_product(kVec, coords(:, iAt)))
      eigVecs(iStart:iEnd, :) = eigVecs(iStart:iEnd, :) * phase
    end do

  end subroutine coordinateIncluding_serial

end module dftbp_math_gauge
