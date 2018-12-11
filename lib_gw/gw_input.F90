module GWInput
#include "allocate.h"  
  use Accuracy
  use SlaterKosterData
  implicit none
  private

  public :: TGWInput, destroy

  type TGWInput
    type(TSlaterKosterData) :: slako
    real(dp), pointer :: xtab(:,:,:) => null()
    real(dp), pointer :: etab(:,:,:) => null()
    real(dp), pointer :: ceri(:,:) => null()
  end type TGWInput


  interface destroy
    module procedure GWInput_destroy
  end interface

contains

  subroutine GWInput_destroy(self)
    type(TGWInput), intent(inout) :: self

    DEALLOCATE_PARR(self%xtab)
    DEALLOCATE_PARR(self%etab)
    DEALLOCATE_PARR(self%ceri)
    
  end subroutine GWInput_destroy
    
end module GWInput
