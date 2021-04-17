!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> input process for GW variables, should replace with slakocont compatible structure
module dftbp_GWInput
  use dftbp_accuracy
  use dftbp_slakocont
  use dftbp_slakoeqgrid
  use dftbp_inputdata
  implicit none
  private

  public :: TGWInput, destroy

  !> slater-koster-like tables of integrals
  type TGWInput
    type(TSlater) :: slako
    real(dp), allocatable :: xtab(:,:,:)
    real(dp), allocatable :: etab(:,:,:)
    real(dp), allocatable :: ceri(:,:)
  end type TGWInput

  !> destroy storage for input data
  interface destroy
    module procedure GWInput_destroy
  end interface

contains

  !> Clean up structure
  subroutine GWInput_destroy(self)

    !> remove structure contents
    type(TGWInput), intent(inout) :: self

    deallocate(self%xtab)
    deallocate(self%etab)
    deallocate(self%ceri)

  end subroutine GWInput_destroy

end module dftbp_GWInput
