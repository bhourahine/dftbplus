!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_math_project
  use dftbp_common_accuracy, only : dp
  use dftbp_math_blasroutines, only : symm
  use dftbp_math_matrixops, only : adjointLowerTriangle
#:if WITH_SCALAPACK
  use dftbp_math_matrixops, only : adjointLowerTriangle_BLACS
#:endif
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: projEmptyReal, projEmptyCplx, projEmptyPauli
#:if WITH_SCALAPACK
  public :: projEmptyReal_BLACS, projEmptyCplx_BLACS, projEmptyPauli_BLACS
#:endif

contains


  subroutine projEmptyReal(P, sSqr, Pc)

    real(dp), intent(in) :: P(:,:,:)

    real(dp), intent(in) :: sSqr(:,:)

    real(dp), intent(out) :: Pc(:,:,:)

    integer :: iSpin, iOrb

    Pc = -P
    do iSpin = 1, size(P, dim=3)
      call adjointLowerTriangle(Pc(:,:,iSpin))
      !call symm(r, 'L', A, x, alpha=-1.0_dp, beta=1.0_dp)
      Pc(:,:,iSpin) = matmul(sSqr,Pc(:,:,iSpin))
      do iOrb = 1, size(Pc, dim=2)
        Pc(iOrb,iOrb,iSpin) = 1.0_dp + Pc(iOrb,iOrb,iSpin)
      end do
    end do
    ! now have (1 - S.Pv)

  end subroutine projEmptyReal


#:if WITH_SCALAPACK

  !> Real projector onto empty states (kernel of the occupied space projector)
  subroutine projEmptyReal_BLACS(grpproc, groupKS, desc, SSqrReal, nFilled, eigvecs, over, species,&
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, Pc)
    type(blacsgrid), intent(in) :: grpproc
    integer, intent(in) :: groupKS(:,:)
    integer, intent(in) :: desc(DLEN_)
    real(dp), intent(inout) :: SSqrReal(:,:)
    integer, intent(in) :: nFilled(:)
    real(dp), intent(in) :: eigvecs(:,:,:)
    real(dp), intent(in) :: over(:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(out) :: Pc(:,:,:)

    integer :: iS, iKS
    integer :: nSpin, nKS, nOrbs, nSparse
    integer :: ii, jj, iGlob, jGlob
    real(dp) :: SPc(size(Pc,dim=1),size(Pc,dim=2),size(Pc,dim=3))

    nSpin = size(eigvecs, dim=3)
    nKS = size(groupKS, dim=2)

    Pc(:,:,:) = 0.0_dp
    SPc(:,:,:) = 0.0_dp

    ! Pc = c . cT
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_psyrk(eigvecs(:,:,iS), desc, Pc(:,:,iS), desc, kk=nFilled(iS))
    end do

    ! fill in other triangle

    ! transpose and add
    SPc = Pc
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_ptran(Pc(:,:,iS),desc,SPc(:,:,iS),desc,beta=1.0_dp)
    end do
    ! correct diagonal, as is double included in transpose above
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      do jj = 1, size(SPc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_), &
            & grpproc%ncol)
        do ii = 1, size(SPc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), grpproc%myrow,desc(RSRC_), &
              & grpproc%nrow)
          if (iGlob == jGlob) then
            SPc(ii,jj,iS) = SPc(ii,jj,iS) * 0.5_dp
          end if
        end do
      end do
    end do

    ! Pc = matmul(SSqr,Pc)
    call unpackhs_parallel_real(grpproc, over, iNeighbor, nNeighbor,&
        & iAtomStart, iPair, img2CentCell, desc, SSqrReal)
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_psymm(SSqrReal, desc, SPc(:,:,iS), desc, Pc(:,:,iS), desc)
    end do

    ! diag(Pc) = diag(Pc) - 1
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      do jj = 1, size(Pc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_), &
            & grpproc%ncol)
        do ii = 1, size(Pc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), grpproc%myrow,desc(RSRC_), &
              & grpproc%nrow)
          if (iGlob == jGlob) then
            Pc(ii,jj,iS) = Pc(ii,jj,iS) - 1.0_dp
          end if
        end do
      end do
    end do

    Pc = -Pc ! now have (1 - S.Pv)

  end subroutine projEmptyReal_BLACS


  !> Complex projector onto conduction band
  subroutine projEmptyCplx_BLACS(grpproc, groupKS, desc, SSqr, nFilled, eigvecs, over, kpoint,&
      & iCellVec, cellVec, species, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, Pc)
    type(blacsgrid), intent(in) :: grpproc
    integer, intent(in) :: groupKS(:,:)
    integer, intent(in) :: desc(DLEN_)
    complex(dp), intent(inout) :: SSqr(:,:)
    integer, intent(in) :: nFilled(:)
    complex(dp), intent(in) :: eigvecs(:,:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kpoint(3)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)
    complex(dp), intent(out) :: Pc(:,:,:)

    integer :: iS, iKS
    integer :: nSpin, nKS, nOrbs, nSparse
    integer :: ii, jj, iGlob, jGlob
    complex(dp) :: SPc(size(Pc,dim=1),size(Pc,dim=2),size(Pc,dim=3))

    nSpin = size(eigvecs, dim=3)
    nKS = size(groupKS, dim=2)

    Pc(:,:,:) = 0.0_dp
    SPc(:,:,:) = 0.0_dp

    ! Pc = c . cT
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_pherk(eigvecs(:,:,iS), desc, Pc(:,:,iS), desc, &
          & kk=nFilled(iS))
    end do

    ! fill in other triangle

    ! Hermitian transpose and add
    SPc = Pc
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_ptranc(Pc(:,:,iS),desc,SPc(:,:,iS),desc,beta=(1.0_dp,0.0_dp))
    end do
    ! correct diagonal, as is double included in transpose above
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      do jj = 1, size(SPc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_), &
            & grpproc%ncol)
        do ii = 1, size(SPc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), grpproc%myrow,desc(RSRC_), &
              & grpproc%nrow)
          if (iGlob == jGlob) then
            SPc(ii,jj,iS) = SPc(ii,jj,iS) * 0.5_dp
          end if
        end do
      end do
    end do

    ! Pc = matmul(SSqr,Pc)
    call unpackhs_parallel_cplx(grpproc, over, kPoint, iNeighbor, nNeighbor, iCellVec, cellVec,&
        & iAtomStart, iPair, img2CentCell, desc, SSqr)
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_phemm(SSqr, desc, SPc(:,:,iS), desc, Pc(:,:,iS), desc)
    end do

    ! diag(Pc) = diag(Pc) - 1
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      do jj = 1, size(Pc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_), &
            & grpproc%ncol)
        do ii = 1, size(Pc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), grpproc%myrow,desc(RSRC_), &
              & grpproc%nrow)
          if (iGlob == jGlob) then
            Pc(ii,jj,iS) = Pc(ii,jj,iS) - 1.0_dp
          end if
        end do
      end do
    end do

    Pc = -Pc ! now have (1 - S.Pv)

  end subroutine projEmptyCplx_BLACS


  !> Pauli matrix projection
  subroutine projEmptyPauli_BLACS(grpproc, groupKS, desc, SSqr, nFilled, eigvecs, over, kpoint,&
      & iCellVec, cellVec, orb, species, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, Pc)
    type(blacsgrid), intent(in) :: grpproc
    integer, intent(in) :: groupKS(:,:)
    integer, intent(in) :: desc(DLEN_)
    complex(dp), intent(inout) :: SSqr(:,:)
    integer, intent(in) :: nFilled(:)
    complex(dp), intent(in) :: eigvecs(:,:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kpoint(3)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)
    complex(dp), intent(out) :: Pc(:,:,:)

    integer :: iS, iKS
    integer :: nSpin, nKS, nOrbs, nSparse
    integer :: ii, jj, iGlob, jGlob
    complex(dp) :: SPc(size(Pc,dim=1),size(Pc,dim=2),size(Pc,dim=3))

    nSpin = size(eigvecs, dim=3)
    nKS = size(groupKS, dim=2)

    Pc(:,:,:) = 0.0_dp
    SPc(:,:,:) = 0.0_dp

    ! Pc = c . cT
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_pherk(eigvecs(:,:,iS), desc, Pc(:,:,iS), desc, &
          & kk=nFilled(iS))
    end do

    ! fill in other triangle

    ! Hermitian transpose and add
    SPc = Pc
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_ptranc(Pc(:,:,iS),desc,SPc(:,:,iS),desc,beta=(1.0_dp,0.0_dp))
    end do
    ! correct diagonal, as is double included in transpose above
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      do jj = 1, size(SPc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_), &
            & grpproc%ncol)
        do ii = 1, size(SPc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), grpproc%myrow,desc(RSRC_), &
              & grpproc%nrow)
          if (iGlob == jGlob) then
            SPc(ii,jj,iS) = SPc(ii,jj,iS) * 0.5_dp
          end if
        end do
      end do
    end do

    ! Pc = matmul(SSqr,Pc)
    call unpacks_parallel_pauli(grpproc, over, kPoint, iNeighbor, &
        & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell, &
        & orb%mOrb, desc, SSqr)
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      call pblasfx_phemm(SSqr, desc, SPc(:,:,iS), desc, Pc(:,:,iS), desc)
    end do

    ! diag(Pc) = diag(Pc) - 1
    do iKS = 1, nKS
      iS = groupKS(2, iKS)
      do jj = 1, size(Pc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_), &
            & grpproc%ncol)
        do ii = 1, size(Pc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), grpproc%myrow,desc(RSRC_), &
              & grpproc%nrow)
          if (iGlob == jGlob) then
            Pc(ii,jj,iS) = Pc(ii,jj,iS) - 1.0_dp
          end if
        end do
      end do
    end do

    Pc = -Pc ! now have (1 - S.Pv)

  end subroutine projEmptyPauli_BLACS

#:endif

end module dftbp_math_project
