!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains BLACS assistance routines
module dftbp_common_blacshelper
  use dftbp_common_accuracy, only : dp
  use dftbp_common_blacsenv, only : TBlacsEnv
  use dftbp_common_memman, only : resizeMat
  use dftbp_extlibs_scalapackfx, only : MB_, NB_, CSRC_, RSRC_, scalafx_indxl2g
  use dftbp_math_sorting, only : heap_sort, unique
  use dftbp_math_bisect, only : bisection
  use dftbp_type_densedescr, only : TDenseDescr
  implicit none

  private
  public :: localAtomsOrbGrid

#:set FLAVOURS = [('Real', 'real'), ('Cmplx', 'complex')]

  interface localAtomsOrbGrid
  #:for SUFFIX, _ in FLAVOURS
    module procedure local${SUFFIX}$
  #:endfor
  end interface localAtomsOrbGrid

contains

#:for SUFFIX, VARTYPE in FLAVOURS

  !> List of atoms associated with local part of an orbital grid matrix
  function local${SUFFIX}$(myBlacs, denseDescr, localMatrix) result (ats)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: denseDescr

    !> Local matrix, used for sizings
    ${VARTYPE}$(dp), intent(in) :: localMatrix(:,:)

    !> List of atoms which have this local matrix as a neighbour
    integer, allocatable :: ats(:)

    integer :: ii, jj, iOrb, jOrb, iAt, jAt, m, n, nAt, mAt, nOrb
    integer, allocatable :: iWork(:)

    m = size(localMatrix, dim=1)
    n = size(localMatrix, dim=2)

    ! maximum number of orbitals (just in case this is indexing an array using spinors)
    ii = size(denseDescr%iAtomStart)
    nOrb = denseDescr%iAtomStart(ii) - 1

    allocate(iWork(max(m,n)), source=0)

    nAt = 0
    do jj = 1, size(localMatrix, dim=2)
      jOrb = scalafx_indxl2g(jj, denseDescr%blacsOrbSqr(NB_), myBlacs%orbitalGrid%mycol,&
          & denseDescr%blacsOrbSqr(CSRC_), myBlacs%orbitalGrid%ncol)
      ! spinor adjustment, if needed:
      jOrb = mod(jOrb-1, nOrb)+1
      call bisection(jAt, denseDescr%iAtomStart, jOrb)
      nAt = nAt + 1
      call resizeMat(iWork, nAt)
      iWork(nAt) = jAt
      do ii = 1, size(localMatrix, dim=1)
        iOrb = scalafx_indxl2g(ii, denseDescr%blacsOrbSqr(MB_), myBlacs%orbitalGrid%myrow,&
            & denseDescr%blacsOrbSqr(RSRC_), myBlacs%orbitalGrid%nrow)
        ! spinor adjustment, if needed:
        iOrb = mod(iOrb-1, nOrb)+1
        call bisection(iAt, denseDescr%iAtomStart, iOrb)
        nAt = nAt + 1
        call resizeMat(iWork, nAt)
        iWork(nAt) = iAt
      end do
    end do

    call heap_sort(iWork(:nAt))
    mAt = unique(iWork(:nAt))
    ats = iWork(:mAt)

  end function local${SUFFIX}$

#:endfor

end module dftbp_common_blacshelper
