!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module containing routines for Sternheimer linear response
module dftbp_derivs_sternheimer
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: sternheimerDynamic, sternheimerStatic

  interface sternheimerDynamic
    module procedure sternheimerDynamic1
    module procedure sternheimerDynamic3
  end interface sternheimerDynamic

contains


!  subroutine sternheimerRef(filling, ei, ci, ham, species)
!    real(dp), intent(in)            :: filling(:,:,:)
!    real(dp), intent(in)            :: ei(:,:,:)
!    real(dp), intent(inout)         :: ci(:,:,:)
!    real(dp), intent(in)            :: ham(:)
!    integer, intent(in)             :: species(:)
!
!    integer  :: ii, iCart, iAtom, iFilled
!    integer  :: nFilled(nSpin)
!    integer  :: iSpin
!    real(dp), allocatable :: hprime(:,:)
!    type(TPotentials) :: dpotential
!    real(dp), allocatable :: HprimeSqr(:,:)
!    real(dp), allocatable :: SSqr(:,:)
!    real(dp), allocatable :: dei(:,:), dci(:,:)
!    real(dp) :: drho(size(over),nSpin)
!    real(dp) :: dqIn(orb%mOrb,nAtom,nSpin)
!    real(dp) :: dqOut(orb%mOrb,nAtom,nSpin)
!    real(dp), allocatable  :: rhoV(:,:)
!    real(dp), allocatable  :: Pc(:,:)
!    real(dp), allocatable  :: PcT(:,:)
!    real(dp), allocatable :: ciTmp(:,:)
!    logical :: tMetallic
!    real(dp) :: polarisability(3,3)
!
!    ALLOCATE_(HprimeSqr,(size(ci,dim=1),size(ci,dim=2)))
!    ALLOCATE_(SSqr,(size(ci,dim=1),size(ci,dim=2)))
!    ALLOCATE_(rhoV,(size(ci,dim=1),size(ci,dim=2)))
!    ALLOCATE_(Pc,(size(ci,dim=1),size(ci,dim=2)))
!    ALLOCATE_(PcT,(size(ci,dim=1),size(ci,dim=2)))
!
!    if (tWriteTagged) then
!      open(fdTagged, file=taggedOut, position="append")
!      call writeTagged(fdTagged, tag_StatPol, .true.)
!    end if
!
!    if (tSCC) call error("Not SCC yet")
!
!    do iSpin = 1, nSpin
!      nFilled(iSpin) = floor(nEl(iSpin)/real(3-nSpin,dp))
!    end do
!
!    call create(dpotential,orb,nAtom,nSpin)
!
!    tMetallic = .false.
!    do iSpin = 1, nSpin
!      if ( any(filling(:nFilled(ispin),1,iSpin) /= real(3-nSpin,dp)) .or. &
!          & any(filling(nFilled(iSpin)+1:,1,iSpin) /= 0.0_dp) ) then
!        tMetallic = .true.
!        call error("Apparently metallic in linear response routine.")
!      end if
!    end do
!
!    ALLOCATE_(hprime,(size(over),nSpin))
!    if (tMetallic) then
!      ALLOCATE_(dei,(size(ci,dim=2),nSpin))
!    else
!      ALLOCATE_(dei,(maxval(nFilled),nSpin))
!    end if
!    ALLOCATE_(dci,(size(ci,dim=1),size(ci,dim=2)))
!    dci = 0.0_dp
!
!    Pc = matmul(ci(:,:nFilled(1),1),transpose(ci(:,:nFilled(1),1)))
!
!    call unpackHS(SSqr,over,neighborList%iNeighbor, &
!        & nNeighbor, iAtomStart,iPair,img2CentCell)
!
!    do ii = 1, size(ci,dim=2)
!      SSqr(ii,ii+1:) = SSqr(ii+1:,ii)
!    end do
!
!    PcT = matmul(Pc,SSqr)
!    Pc = matmul(SSqr,Pc)
!
!    rhoV = matmul(SSqr,PcT)
!
!    do ii = 1, size(ci,dim=2)
!      Pc(ii,ii) = Pc(ii,ii) - 1.0_dp
!      PcT(ii,ii) = PcT(ii,ii) - 1.0_dp
!    end do
!
!    Pc = - Pc ! (1 - S.Pv)
!    PcT = - PcT ! (1 - Pv.S)
!
!    ALLOCATE_(ciTmp,(size(ci,dim=1),maxval(nFilled)))
!
!    do iCart = 1, 3
!
!      dqIn = 0.0_dp
!
!      dpotential%extAtom = 0.0_dp
!      do iAtom = 1, nAtom
!        dpotential%extAtom(iAtom,1) = coord0(iCart,iAtom)
!      end do
!
!      dpotential%extBlock = 0.0_dp
!      call copy_shiftTo(dpotential%extBlock,dpotential%extAtom, &
!          & orb,species)
!
!      hprime = 0.0_dp
!
!      call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
!          & species,orb,iPair,nAtom,img2CentCell,dpotential%extBlock)
!      dei = 0.0_dp
!      HprimeSqr = 0.0_dp
!      call unpackHS(HprimeSqr,hprime(:,1),neighborList%iNeighbor, &
!          & nNeighbor, iAtomStart,iPair,img2CentCell)
!
!      do ii = 1, size(ci,dim=2)
!        HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
!      end do
!
!      ciTmp(:,:nFilled(1)) = matmul(HprimeSqr,ci(:,:nFilled(1),1))
!      ciTmp(:,:nFilled(1)) = -matmul(Pc,ciTmp(:,:nFilled(1)))
!
!      do iFilled = 1, nFilled(1)
!
!        HprimeSqr = 0.0_dp
!        call unpackHS(HprimeSqr,ham,neighborList%iNeighbor, &
!            & nNeighbor, iAtomStart,iPair,img2CentCell)
!
!        SSqr = 0.0_dp
!        call unpackHS(SSqr,over,neighborList%iNeighbor, &
!            & nNeighbor, iAtomStart,iPair,img2CentCell)
!
!        HprimeSqr = HprimeSqr - ei(iFilled,1,1)*SSqr
!
!        do ii = 1, size(ci,dim=2)
!          HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
!        end do
!
!        dci(:,iFilled) = 1.0_dp
!        call conjgrad(HprimeSqr,ciTmp(:,iFilled),dci(:,iFilled))
!
!      end do
!
!      dci(:,:nFilled(1)) = matmul(PcT,dci(:,:nFilled(1)))
!
!      dci = ( matmul(dci(:,:nFilled(1)), &
!          & transpose(ci(:,:nFilled(1),1))) + &
!          & matmul(ci(:,:nFilled(1),1), &
!          & transpose(dci(:,:nFilled(1)))) )
!
!      drho(:,1) = 0.0_dp
!      call packHS(drho(:,1), dci, neighborlist%iNeighbor, &
!          & nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
!
!      drho = real(3-nSpin,dp) * drho
!
!      dqOut = 0.0_dp
!      do iSpin = 1, nSpin
!        call mulliken(dqOut(:,:,iSpin), over, drho(:,iSpin), orb, &
!            & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
!      end do
!
!      if (iCart==1) then
!        write(*,*)'Polarisability'
!      end if
!      do ii = 1, 3
!        write(*,"(E20.12)",advance ='no')-sum(sum(dqOut(:,:,1),dim=1)*coord0(ii&
!            &,:))
!      end do
!      write(*,*)
!
!    end do
!
!    call destroy(dpotential)
!    DEALLOCATE_(ciTmp)
!    DEALLOCATE_(dci)
!    DEALLOCATE_(dei)
!    DEALLOCATE_(hprime)
!    DEALLOCATE_(HprimeSqr)
!    DEALLOCATE_(SSqr)
!    DEALLOCATE_(rhoV)
!    DEALLOCATE_(Pc)
!    DEALLOCATE_(PcT)
!
!
!  end subroutine sternheimerRef

  subroutine sternheimerStatic(filling, ei, ci, ham, species, tDense, tol)
    real(dp), intent(in)            :: filling(:,:,:)
    real(dp), intent(in)            :: ei(:,:,:)
    real(dp), intent(inout)         :: ci(:,:,:)
    real(dp), intent(in)            :: ham(:,:)
    integer, intent(in)             :: species(:)
    logical, intent(in)             :: tDense
    real(dp), intent(in)            :: tol

    integer  :: ii, iCart, iAtom
    integer  :: nFilled(nSpin)
    integer  :: iSpin, iSpin2, iSh1, iSp1
    real(dp), allocatable :: hprime(:,:)
    type(TPotentials) :: dpotential
    real(dp), allocatable :: HprimeSqr(:,:)
    real(dp), allocatable :: SSqr(:,:)
    real(dp), allocatable :: dei(:,:), dci(:,:,:)
    real(dp) :: drho(size(over),nSpin)
    real(dp) :: dqIn(orb%mOrb,nAtom,nSpin)
    real(dp) :: dqOut(orb%mOrb,nAtom,nSpin)
    real(dp), allocatable :: Pc(:,:,:), arrayTmp(:,:), arrayTmp2(:,:)
    real(dp), allocatable :: PcT(:,:,:)
    real(dp), allocatable :: ciTmp(:,:,:)

    real(dp) :: sccErrorQ, polarisability(3,3)
    integer  :: iSCCIter
    logical  :: tStopSCC, tConverged, tMetallic
    character, parameter :: directLabel(3) = ['x','y','z']

    if (.not.tDense) then
      write(*,*)'*****************************************************'
      write(*,*)'Sparse implementation (partially)'
      write(*,*)'*****************************************************'
    end if

    if (tWriteTagged) then
      open(fdTagged, file=taggedOut, position="append")
      call writeTagged(fdTagged, tag_StatPol, .true.)
    end if

    ALLOCATE_(HprimeSqr,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(SSqr,(size(ci,dim=1),size(ci,dim=2)))

    ALLOCATE_(Pc,(size(ci,dim=1),size(ci,dim=2),nSpin))
    ALLOCATE_(PcT,(size(ci,dim=1),size(ci,dim=2),nSpin))

    do iSpin = 1, nSpin
      nFilled(iSpin) = floor(nEl(iSpin)/real(3-nSpin,dp))
    end do

    call create(dpotential,orb,nAtom,nSpin)

    tMetallic = .false.
    do iSpin = 1, nSpin
      if ( any(filling(:nFilled(ispin),1,iSpin) /= real(3-nSpin,dp)) .or. &
          & any(filling(nFilled(iSpin)+1:,1,iSpin) /= 0.0_dp) ) then
        tMetallic = .true.
        call error("Apparently metallic in static Sternheimer routine.")
      end if
    end do

    ALLOCATE_(hprime,(size(over),nSpin))
    ALLOCATE_(dei,(maxval(nFilled),nSpin))
    ALLOCATE_(dci,(size(ci,dim=1),size(ci,dim=2),nSpin))

    write(*,*)'Forming Pc'

    ! Pc = matmul(ci(:,:nFilled(1),1),transpose(ci(:,:nFilled(1),1)))
    Pc = 0.0_dp
    do iSpin = 1, nSpin
      call herk(Pc(:,:,iSpin),ci(:,:nFilled(1),iSpin))
      do ii = 1, size(ci,dim=2)
        Pc(ii,ii+1:,iSpin) = Pc(ii+1:,ii,iSpin)
      end do
    end do

    SSqr = 0.0_dp
    call unpackHS(SSqr,over,neighborList%iNeighbor, &
        & nNeighbor, iAtomStart,iPair,img2CentCell)
    do ii = 1, size(ci,dim=2)
      SSqr(ii,ii+1:) = SSqr(ii+1:,ii)
    end do

    write(*,*)'Forming Pc and PcT'
    do iSpin = 1, nSpin
      call gemm(Pct(:,:,iSpin),Pc(:,:,iSpin),SSqr,transA='T')
    end do
    !PcT = matmul(Pc,SSqr)
    ALLOCATE_(arrayTmp,(size(ci,dim=1),size(ci,dim=2)))
    !Pc = matmul(SSqr,Pc)
    do iSpin = 1, nSpin
      call gemm(arrayTmp,SSqr,Pc(:,:,iSpin),transA='T')
      Pc(:,:,iSpin) = arrayTmp
    end do
    DEALLOCATE_(arrayTmp)

    do iSpin = 1, nSpin
      do ii = 1, size(ci,dim=2)
        Pc(ii,ii,iSpin) = Pc(ii,ii,iSpin) - 1.0_dp
        PcT(ii,ii,iSpin) = PcT(ii,ii,iSpin) - 1.0_dp
      end do
    end do

    Pc = - Pc ! (1 - S.Pv)
    PcT = - PcT ! (1 - Pv.S)

    write(*,*)'Now have projectors'

    ALLOCATE_(ciTmp,(size(ci,dim=1),maxval(nFilled),nSpin))

    do iCart = 1, 3
      write(*,*)'Direction ',directLabel(iCart)
      dci = 1.0_dp
      dqIn = 0.0_dp

      dpotential%extAtom = 0.0_dp
      do iAtom = 1, nAtom
        dpotential%extAtom(iAtom,1) = coord0(iCart,iAtom)
      end do

      dpotential%extBlock = 0.0_dp
      call copy_shiftTo(dpotential%extBlock,dpotential%extAtom, &
          & orb,species)
      !! (Re)Initialize mixer
      if (tSCC) then
        call reset(pChrgMixer,nMixElements) !nIneqOrb)
        chargePerShell = 0.0_dp
        qBlockIn = 0.0_dp
        qInpRed = 0.0_dp
      end if
      iSCCIter = 1
      tStopSCC = .false.
      lpSCC: do while (iSCCiter <= nSCCIter)

        dpotential%intAtom = 0.0_dp
        dpotential%intBlock = 0.0_dp
        dpotential%intShell = 0.0_dp
        dpotential%orbitalBlock = 0.0_dp

        dpotential%orbitalBlock = dpotential%extBlock
        if (tSCC .and. iSCCiter>1) then
          call updateCharges_SCC(dqIn, orb, species, &
              &neighborList%iNeighbor, img2CentCell)
          call getShiftPerAtom(dpotential%intAtom)
          call getShiftPerL(dpotential%intShell)

          call total_shift(dpotential%intShell,dpotential%intAtom, &
              & dpotential%intShell,orb,species)

          if (tSpin) then
            call addSpinShift(dpotential%intShell,chargePerShell,species,orb,W)
          end if
          if (tDFTBU) then
            call shift_DFTBU(dpotential%orbitalBlock,qBlockIn,specie,orb, &
                & 2, UJ, nUJ, niUJ, iUJ)
          end if
        end if

        call total_shift(dpotential%intBlock,dpotential%intShell, &
            & dpotential%orbitalBlock,orb,species)

        hprime = 0.0_dp

        call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
            & species,orb,iPair,nAtom,img2CentCell,dpotential%intBlock)

        if (nSpin == 2) then
          hprime(:,:) = 2.0_dp *  hprime(:,:)
        end if
        call qm2ud(hprime)

        dei = 0.0_dp

        do iSpin = 1, nSpin

          HprimeSqr = 0.0_dp
          call unpackHS(HprimeSqr,hprime(:,iSpin),neighborList%iNeighbor, &
              & nNeighbor, iAtomStart,iPair,img2CentCell)

          do ii = 1, size(ci,dim=2)
            HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
          end do

          if (tStoreEigvecs) then
            iSpin2 = 1
            call get(storeEigvecsReal(iSpin), &
                & ci(:,:,iSpin2))
          else
            iSpin2 = iSpin
          end if

          dei(1:nFilled(iSpin),iSpin) = sum(ci(:,1:nFilled(iSpin),iSpin2) &
              & *matmul(HprimeSqr,ci(:,1:nFilled(iSpin),iSpin2)),dim=1)

          write(450,*)'Iter', iSCCiter
          write(450,*)'dEi (eV)'
          do ii = 1, nFilled(iSpin)
            write(450,*)dei(ii,iSpin)*Hartree__eV
          end do
          write(450,*)

          call gemm(ciTmp(:,:nFilled(iSpin),iSpin),HprimeSqr, &
              & ci(:,:nFilled(iSpin),iSpin2),transA='T')
          !ciTmp(:,:nFilled(1)) = matmul(HprimeSqr,ci(:,:nFilled(1),1))

          ALLOCATE_(arrayTmp,(size(ciTmp,dim=1),nFilled(iSpin)))

          arrayTmp = 0.0_dp

          !ciTmp(:,:nFilled(1)) = -matmul(Pc,ciTmp(:,:nFilled(1)))
          call gemm(arrayTmp,Pc(:,:,iSpin),ciTmp(:,:nFilled(iSpin),iSpin))
          !arrayTmp  = matmul(Pc,ciTmp(:,:nFilled(1)))
          ciTmp(:,:nFilled(iSpin),iSpin) = -arrayTmp(:,:nFilled(iSpin))
          DEALLOCATE_(arrayTmp)

          if (tDense) then

            HprimeSqr = 0.0_dp
            call unpackHS(HprimeSqr,ham(:,iSpin),neighborList%iNeighbor, &
                & nNeighbor, iAtomStart,iPair,img2CentCell)

            do ii = 1, size(ci,dim=2)
              HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
            end do

            !write(*,*)'entering dense CG routine'
            call conjgrad_contract(HprimeSqr,ei(:nFilled(iSpin),1,iSpin),SSqr, &
                &  ciTmp(:,:nFilled(iSpin),iSpin),dci(:,:nFilled(iSpin),iSpin),&
                & tol=tol)
            !write(*,*)'Leaving CG'
          else

            !write(*,*)'entering sparse CG routine'
            call conjgrad_contract(ham(:,iSpin),ei(:nFilled(iSpin),1,iSpin),over, &
                & ciTmp(:,:nFilled(iSpin),iSpin),dci(:,:nFilled(iSpin),iSpin), &
                & neighborList%iNeighbor, nNeighbor, img2CentCell,iPair, &
                & iAtomStart, orb%mOrb,tol=tol)
            !write(*,*)'Leaving CG'

          end if

          ALLOCATE_(arrayTmp,(size(dci,dim=1),nFilled(iSpin)))

          call gemm(arrayTmp,PcT(:,:,iSpin),dci(:,:nFilled(iSpin),iSpin))

          dci(:,:nFilled(iSpin),iSpin) = arrayTmp

          !dci(:,:nFilled(1)) = matmul(PcT,dci(:,:nFilled(1)))

          !dci = ( matmul(dci(:,:nFilled(1)), &
          !    & transpose(ci(:,:nFilled(1),1))) + &
          !    & matmul(ci(:,:nFilled(1),1), &
          !    & transpose(dci(:,:nFilled(1)))) )

          ALLOCATE_(arrayTmp2,(size(dci,dim=1),size(dci,dim=1)))

          arrayTmp = ci(:,:nFilled(iSpin),iSpin)

          call gemm(arrayTmp2,dci(:,:nFilled(iSpin),iSpin),arrayTmp,transB='T')

          arrayTmp = dci(:,:nFilled(iSpin),iSpin)

          call gemm(arrayTmp2,ci(:,:nFilled(iSpin),iSpin),arrayTmp, &
              & transB='T',beta=1.0_dp)
          dci(:,:,iSpin) = arrayTmp2
          DEALLOCATE_(arrayTmp)
          DEALLOCATE_(arrayTmp2)

          drho(:,iSpin) = 0.0_dp
          call packHS(drho(:,iSpin), dci(:,:,iSpin), neighborlist%iNeighbor, &
              & nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)

        end do

        if (nSpin == 1) then
          drho = 2.0_dp * drho
        end if

        call ud2qm(drho)

        dqOut = 0.0_dp
        do iSpin = 1, nSpin
          call mulliken(dqOut(:,:,iSpin), over, drho(:,iSpin), orb, &
              & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
          if (tDFTBU) then
            qBlockOut(:,:,:,iSpin) = 0.0_dp
            call mulliken(qBlockOut(:,:,:,iSpin), over, drho(:,iSpin), &
                &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
          end if
        end do

        if (tSCC) then
          write(lcTmp,"('chargeDeriv',i1,'.bin')")iCart
          call writeQToFile(dqOut,trim(lcTmp),orb)
        end if


        if (tSCC) then
          qOutRed = 0.0_dp
          call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
              & qOutRed(1:nIneqOrb))
          if (tDFTBU) then
            call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
                & qOutRed )
          end if

          qDiffRed(:) = qOutRed(:) - qInpRed(:)
          sccErrorQ = maxval(abs(qDiffRed))

          write(*,*)iSCCIter,sccErrorQ
          tConverged = (sccErrorQ < sccTol)

          if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
            if (iSCCIter == 1) then
              dqIn(:,:,:) = dqOut(:,:,:)
              qInpRed(:) = qOutRed(:)
              if (tDFTBU) then
                qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
              end if
            else
              call mix(pChrgMixer, qInpRed, qDiffRed)

              call OrbitalEquiv_expand(qInpRed(1:nIneqOrb), iEqOrbitals, &
                  & orb, dqIn)
              if (tDFTBU) then
                qBlockIn = 0.0_dp
                call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
                    & qBlockIn, specie0, nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
              end if
            end if
          end if
        else
          tConverged = .true.
        end if

        write(7001,"(1x,A,I4,A,E12.6)")directLabel(iCart), iSccIter, ' Err ',sccErrorQ
        do ii = 1, 3
          write(7001,"(E20.12)",advance ='no')-sum(sum(dqOut(:,:,1),dim=1)*coord0(ii,:))
        end do
        write(7001,*)
        call flush(7001)


        if (tConverged) then
          exit lpSCC
        end if

        if (tSpin) then
          chargePerShell(:,:,:) = 0.0_dp
          do iAtom = 1, nAtom
            iSp1 = species(iAtom)
            do iSh1 = 1, orb%nShell(iSp1)
              chargePerShell(iSh1,iAtom,1:nSpin) = &
                  & chargePerShell(iSh1,iAtom,1:nSpin) + &
                  & sum(dqIn(orb%posShell(iSh1,iSp1): &
                  & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
            end do
          end do

        end if

        iSCCIter = iSCCIter +1


      end do lpSCC

      !if (iCart==1) then
      write(*,*)'Polarisability ',directLabel(iCart)
      !end if
      do ii = 1, 3
        write(*,"(E20.12)",advance ='no')-sum(sum(dqOut(:,:,1),dim=1)*coord0(ii&
            &,:))
        polarisability(ii,iCart) = -sum(sum(dqOut(:,:,1),dim=1)*coord0(ii,:))
      end do
      write(*,*)

      !do iAtom = 1, nAtom
      !  write(*,*)dqOut(:,iAtom,1)
      !end do
      !write(*,*)

    end do

    call destroy(dpotential)
    DEALLOCATE_(ciTmp)
    DEALLOCATE_(dci)
    DEALLOCATE_(dei)
    DEALLOCATE_(hprime)

    DEALLOCATE_(HprimeSqr)
    DEALLOCATE_(SSqr)
    DEALLOCATE_(Pc)
    DEALLOCATE_(PcT)

    if (tWriteTagged) then
      call writeTagged(fdTagged, tag_StatEPol, polarisability)
      close(fdTagged)
    end if

  end subroutine sternheimerStatic

  ! Dense dynamical Sternheimer
  subroutine sternheimerDynamic1(filling, ei, ci, ham, species, omega)
    real(dp), intent(in)            :: filling(:,:,:)
    real(dp), intent(in)            :: ei(:,:,:)
    real(dp), intent(inout)         :: ci(:,:,:)
    real(dp), intent(in)            :: ham(:)
    integer, intent(in)             :: species(:)
    real(dp), intent(in)            :: omega

    integer  :: ii, jj, iCart, iAtom, iFilled
    integer  :: nFilled(nSpin)
    integer  :: iSpin
    real(dp), allocatable :: hprime(:,:)
    type(TPotentials) :: dpotential
    real(dp), allocatable :: HprimeSqr(:,:)
    real(dp), allocatable :: SSqr(:,:)
    real(dp), allocatable :: HprimeTmp(:,:)
    real(dp), allocatable :: dei(:,:), dci(:,:)
    real(dp) :: drho(size(over),nSpin)
    real(dp) :: dqIn(orb%mOrb,nAtom,nSpin)
    real(dp) :: dqOut(orb%mOrb,nAtom,nSpin)
    real(dp), allocatable  :: rhoV(:,:)
    real(dp), allocatable  :: Pc(:,:)
    real(dp), allocatable  :: PcT(:,:)
    real(dp), allocatable :: ciTmp(:,:)
    logical :: tMetallic

    write(*,*)'Dynamic Sternheimer'

    ALLOCATE_(HprimeSqr,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(SSqr,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(rhoV,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(Pc,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(PcT,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(HprimeTmp,(size(ci,dim=1),size(ci,dim=2)))

    if (tWriteTagged) then
      open(fdTagged, file=taggedOut, position="append")
      call writeTagged(fdTagged, tag_StatPol, .true.)
    end if

    if (tSCC) call error("Not SCC yet")

    do iSpin = 1, nSpin
      nFilled(iSpin) = floor(nEl(iSpin)/real(3-nSpin,dp))
    end do

    call create(dpotential,orb,nAtom,nSpin)

    tMetallic = .false.
    do iSpin = 1, nSpin
      if ( any(filling(:nFilled(ispin),1,iSpin) /= real(3-nSpin,dp)) .or. &
          & any(filling(nFilled(iSpin)+1:,1,iSpin) /= 0.0_dp) ) then
        tMetallic = .true.
        call error("Apparently metallic in linear response routine.")
      end if
    end do

    ALLOCATE_(hprime,(size(over),nSpin))
    ALLOCATE_(dei,(maxval(nFilled),nSpin))
    ALLOCATE_(dci,(size(ci,dim=1),size(ci,dim=2)))
    dci = 0.0_dp

    Pc = matmul(ci(:,:nFilled(1),1),transpose(ci(:,:nFilled(1),1)))

    call unpackHS(SSqr,over,neighborList%iNeighbor, &
        & nNeighbor, iAtomStart,iPair,img2CentCell)

    do ii = 1, size(ci,dim=2)
      SSqr(ii,ii+1:) = SSqr(ii+1:,ii)
    end do

    PcT = matmul(Pc,SSqr)
    Pc = matmul(SSqr,Pc)

    rhoV = matmul(SSqr,PcT)

    do ii = 1, size(ci,dim=2)
      Pc(ii,ii) = Pc(ii,ii) - 1.0_dp
      PcT(ii,ii) = PcT(ii,ii) - 1.0_dp
    end do

    Pc = - Pc ! (1 - S.Pv)
    PcT = - PcT ! (1 - Pv.S)

    ALLOCATE_(ciTmp,(size(ci,dim=1),maxval(nFilled)))

    do iCart = 1, 3

      dqIn = 0.0_dp

      dpotential%extAtom = 0.0_dp
      do iAtom = 1, nAtom
        dpotential%extAtom(iAtom,1) = coord0(iCart,iAtom)
      end do

      dpotential%extBlock = 0.0_dp
      call copy_shiftTo(dpotential%extBlock,dpotential%extAtom, &
          & orb,species)

      hprime = 0.0_dp

      call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
          & species,orb,iPair,nAtom,img2CentCell,dpotential%extBlock)
      dei = 0.0_dp
      HprimeSqr = 0.0_dp
      call unpackHS(HprimeSqr,hprime(:,1),neighborList%iNeighbor, &
          & nNeighbor, iAtomStart,iPair,img2CentCell)

      do ii = 1, size(ci,dim=2)
        HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
      end do

      ciTmp(:,:nFilled(1)) = matmul(HprimeSqr,ci(:,:nFilled(1),1))
      ciTmp(:,:nFilled(1)) = -matmul(Pc,ciTmp(:,:nFilled(1)))

      SSqr = 0.0_dp
      call unpackHS(SSqr,over,neighborList%iNeighbor, &
          & nNeighbor, iAtomStart,iPair,img2CentCell)
      do ii = 1, size(ci,dim=2)
        SSqr(ii,ii+1:) =  SSqr(ii+1:,ii)
      end do

      HprimeSqr = 0.0_dp
      call unpackHS(HprimeSqr,ham,neighborList%iNeighbor, &
          & nNeighbor, iAtomStart,iPair,img2CentCell)

      do ii = 1, size(ci,dim=2)
        HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
      end do

      drho(:,1) = 0.0_dp

      do jj = -1, 1, 2

        do iFilled = 1, nFilled(1)
          dci(:,iFilled) = 1.0_dp
          dci(:,iFilled) = matmul(PcT,dci(:,iFilled))

          HprimeTmp = HprimeSqr + (jj*omega - ei(iFilled,1,1)) *SSqr
          !call bicgstab( HprimeTmp, ciTmp(:,iFilled), dci(:,iFilled), &
          !    & epsilon(1.0_dp))
          call conjgrad(HprimeTmp, ciTmp(:,iFilled), dci(:,iFilled),&
              & epsilon(1.0_dp))
          !dci(:,iFilled) = matmul(PcT,dci(:,iFilled))
          !write(*,*)jj,iFilled, &
          !    &sum((ciTmp(1,iFilled)-matmul(HprimeTmp,dci(:&
          !    &,iFilled)))*(ciTmp(1,iFilled)-matmul(HprimeTmp,dci(:&
          !    &,iFilled))))

        end do

        dci(:,:nFilled(1)) = matmul(PcT,dci(:,:nFilled(1)))

        dci = matmul(ci(:,:nFilled(1),1), &
            & transpose((dci(:,:nFilled(1)))))
        dci = 0.5_dp * (dci + transpose((dci)))

        call packHS(drho(:,1), real(dci,dp), neighborlist%iNeighbor, &
            & nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)

      end do

      drho = real(3-nSpin,dp) * drho

      dqOut = 0.0_dp
      do iSpin = 1, nSpin
        call mulliken(dqOut(:,:,iSpin), over, drho(:,iSpin), orb, &
            & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
      end do

      if (iCart==1) then
        write(*,*)'Polarisability'
      end if
      do ii = 1, 3
        write(*,"(E20.12)",advance ='no')-sum(sum(dqOut(:,:,1),dim=1)*coord0(ii&
            &,:))
      end do
      write(*,*)

    end do

    call destroy(dpotential)

    DEALLOCATE_(HPrimeTmp)
    DEALLOCATE_(ciTmp)
    DEALLOCATE_(dci)
    DEALLOCATE_(dei)
    DEALLOCATE_(hprime)
    DEALLOCATE_(HprimeSqr)
    DEALLOCATE_(SSqr)
    DEALLOCATE_(rhoV)
    DEALLOCATE_(Pc)
    DEALLOCATE_(PcT)

  end subroutine sternheimerDynamic1

  ! sparse partly
  subroutine sternheimerDynamic2(filling, ei, ci, ham, species, &
      & omega, tol)
    real(dp), intent(in)            :: filling(:,:,:)
    real(dp), intent(in)            :: ei(:,:,:)
    real(dp), intent(inout)         :: ci(:,:,:)
    real(dp), intent(in)            :: ham(:,:)
    integer, intent(in)             :: species(:)
    real(dp), intent(in)            :: omega(:)
    real(dp), intent(in)            :: tol

    integer  :: ii, jj, kk, iCart, iAtom
    integer  :: nFilled(nSpin)
    integer  :: iSpin, iSpin2, iSh1, iSp1
    real(dp), allocatable :: hprime(:,:)
    type(TPotentials) :: dpotential
    real(dp), allocatable :: HprimeSqr(:,:)
    real(dp), allocatable :: SSqr(:,:)
    real(dp), allocatable :: dei(:,:), dci(:,:,:,:)
    real(dp) :: drho(size(over),nSpin)
    real(dp) :: dqIn(orb%mOrb,nAtom,nSpin)
    real(dp) :: dqOut(orb%mOrb,nAtom,nSpin)
    real(dp), allocatable :: Pc(:,:,:), arrayTmp(:,:), arrayTmp2(:,:)
    real(dp), allocatable :: PcT(:,:,:)
    real(dp), allocatable :: ciTmp(:,:,:), eiTmp(:)

    real(dp) :: sccErrorQ, polarisability(3,3)
    integer  :: iSCCIter, iOmega, nOmega
    logical  :: tStopSCC, tConverged, tMetallic

    nOmega = size(omega)

    if (tWriteTagged) then
      open(fdTagged, file=taggedOut, position="append")
      call writeTagged(fdTagged, tag_StatPol, .true.)
    end if

    ALLOCATE_(HprimeSqr,(size(ci,dim=1),size(ci,dim=2)))
    ALLOCATE_(SSqr,(size(ci,dim=1),size(ci,dim=2)))

    ALLOCATE_(Pc,(size(ci,dim=1),size(ci,dim=2),nSpin))
    ALLOCATE_(PcT,(size(ci,dim=1),size(ci,dim=2),nSpin))

    do iSpin = 1, nSpin
      nFilled(iSpin) = floor(nEl(iSpin)/real(3-nSpin,dp))
    end do

    call create(dpotential,orb,nAtom,nSpin)

    tMetallic = .false.
    do iSpin = 1, nSpin
      if ( any(filling(:nFilled(ispin),1,iSpin) /= real(3-nSpin,dp)) .or. &
          & any(filling(nFilled(iSpin)+1:,1,iSpin) /= 0.0_dp) ) then
        tMetallic = .true.
        call error("Apparently metallic in static Sternheimer routine.")
      end if
    end do

    ALLOCATE_(hprime,(size(over),nSpin))
    ALLOCATE_(dei,(maxval(nFilled),nSpin))
    ALLOCATE_(dci,(size(ci,dim=1),size(ci,dim=2),2,nSpin))

    write(*,*)'Forming Pc'

    ! Pc = matmul(ci(:,:nFilled(1),1),transpose(ci(:,:nFilled(1),1)))
    Pc = 0.0_dp
    do iSpin = 1, nSpin
      call herk(Pc(:,:,iSpin),ci(:,:nFilled(1),iSpin))
      do ii = 1, size(ci,dim=2)
        Pc(ii,ii+1:,iSpin) = Pc(ii+1:,ii,iSpin)
      end do
    end do

    SSqr = 0.0_dp
    call unpackHS(SSqr,over,neighborList%iNeighbor, &
        & nNeighbor, iAtomStart,iPair,img2CentCell)
    do ii = 1, size(ci,dim=2)
      SSqr(ii,ii+1:) = SSqr(ii+1:,ii)
    end do

    write(*,*)'Forming Pc and PcT'
    do iSpin = 1, nSpin
      call gemm(Pct(:,:,iSpin),Pc(:,:,iSpin),SSqr,transA='T')
    end do
    !PcT = matmul(Pc,SSqr)
    ALLOCATE_(arrayTmp,(size(ci,dim=1),size(ci,dim=2)))
    !Pc = matmul(SSqr,Pc)
    do iSpin = 1, nSpin
      call gemm(arrayTmp,SSqr,Pc(:,:,iSpin),transA='T')
      Pc(:,:,iSpin) = arrayTmp
    end do
    DEALLOCATE_(arrayTmp)

    do iSpin = 1, nSpin
      do ii = 1, size(ci,dim=2)
        Pc(ii,ii,iSpin) = Pc(ii,ii,iSpin) - 1.0_dp
        PcT(ii,ii,iSpin) = PcT(ii,ii,iSpin) - 1.0_dp
      end do
    end do

    Pc = - Pc ! (1 - S.Pv)
    PcT = - PcT ! (1 - Pv.S)

    write(*,*)'Now have projectors'

    ALLOCATE_(ciTmp,(size(ci,dim=1),maxval(nFilled),nSpin))
    ALLOCATE_(eiTmp,(maxval(nFilled)))
    eiTmp = 0.0_dp

    do iCart = 1, 3

      write(*,*)'Polarisability along ',direction(iCart)

      dci = 1.0_dp
      dqIn = 0.0_dp

      dpotential%extAtom = 0.0_dp
      do iAtom = 1, nAtom
        dpotential%extAtom(iAtom,1) = coord0(iCart,iAtom)
      end do

      dpotential%extBlock = 0.0_dp
      call copy_shiftTo(dpotential%extBlock,dpotential%extAtom, &
          & orb,species)

      do iOmega = 1, nOmega

!        dpotential%intAtom = 0.0_dp
!        dpotential%intShell = 0.0_dp
!        dpotential%intBlock = 0.0_dp
!        dpotential%extAtom = 0.0_dp
!        dpotential%extShell = 0.0_dp
!        dpotential%orbitalBlock = 0.0_dp
!        dpotential%iorbitalBlock = 0.0_dp

!        chargePerShell = 0.0_dp
!        qBlockIn = 0.0_dp
!        qInpRed = 0.0_dp
!        dqIn = 0.0_dp
        dqOut = 0.0_dp

        write(*,*)'At',omega(iOmega),'au',omega(iOmega)*Hartree__eV,'eV'
 !       dci = 1.0_dp


        !! (Re)Initialize mixer
        if (tSCC) then
          call reset(pChrgMixer,nMixElements) !nIneqOrb)
          if (iOmega == 1) then
            chargePerShell = 0.0_dp
            qBlockIn = 0.0_dp
            qInpRed = 0.0_dp
          end if
        end if
        iSCCIter = 1
        tStopSCC = .false.

        lpSCC: do while (iSCCiter <= nSCCIter)

          dpotential%intAtom = 0.0_dp
          dpotential%intBlock = 0.0_dp
          dpotential%intShell = 0.0_dp
          dpotential%orbitalBlock = 0.0_dp
          dpotential%extAtom = 0.0_dp

          dpotential%orbitalBlock = dpotential%extBlock
          if (tSCC .and. (iSCCiter>1 .or. iOmega > 1)) then
            call updateCharges_SCC(dqIn, orb, species, &
                &neighborList%iNeighbor, img2CentCell)
            call getShiftPerAtom(dpotential%intAtom)
            call getShiftPerL(dpotential%intShell)

            call total_shift(dpotential%intShell,dpotential%intAtom, &
                & dpotential%intShell,orb,species)

            if (tSpin) then
              call addSpinShift(dpotential%intShell,chargePerShell,species,orb&
                  &,W)
            end if
            if (tDFTBU) then
              call shift_DFTBU(dpotential%orbitalBlock,qBlockIn,specie,orb, &
                  & 2, UJ, nUJ, niUJ, iUJ)
            end if
          end if

          call total_shift(dpotential%intBlock,dpotential%intShell, &
              & dpotential%orbitalBlock,orb,species)

          hprime = 0.0_dp

          call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
              & species,orb,iPair,nAtom,img2CentCell,dpotential%intBlock)

          if (nSpin == 2) then
            hprime(:,:) = 2.0_dp *  hprime(:,:)
          end if
          call qm2ud(hprime)

          dei = 0.0_dp
          drho = 0.0_dp

          do iSpin = 1, nSpin

            HprimeSqr = 0.0_dp
            call unpackHS(HprimeSqr,hprime(:,iSpin),neighborList%iNeighbor, &
                & nNeighbor, iAtomStart,iPair,img2CentCell)

            do ii = 1, size(ci,dim=2)
              HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
            end do

            if (tStoreEigvecs) then
              iSpin2 = 1
              call get(storeEigvecsReal(iSpin), ci(:,:,iSpin2))
            else
              iSpin2 = iSpin
            end if

            dei(1:nFilled(iSpin),iSpin) = sum(ci(:,1:nFilled(iSpin),iSpin2) &
                & *matmul(HprimeSqr,ci(:,1:nFilled(iSpin),iSpin2)),dim=1)

            call gemm(ciTmp(:,:nFilled(iSpin),iSpin),HprimeSqr, &
                & ci(:,:nFilled(iSpin),iSpin2),transA='T')
            !ciTmp(:,:nFilled(1)) = matmul(HprimeSqr,ci(:,:nFilled(1),1))

            ALLOCATE_(arrayTmp,(size(ciTmp,dim=1),nFilled(iSpin)))

            arrayTmp = 0.0_dp

            !ciTmp(:,:nFilled(1)) = -matmul(Pc,ciTmp(:,:nFilled(1)))
            call gemm(arrayTmp,Pc(:,:,iSpin),ciTmp(:,:nFilled(iSpin),iSpin))
            !arrayTmp  = matmul(Pc,ciTmp(:,:nFilled(1)))
            ciTmp(:,:nFilled(iSpin),iSpin) = -arrayTmp(:,:nFilled(iSpin))
            DEALLOCATE_(arrayTmp)


            do jj = 1, 2
              kk = 2*jj-3


              write(*,*)'entering sparse CG routine'
              eiTmp(:nFilled(iSpin)) = &
                  & ei(:nFilled(iSpin),1,iSpin)+kk*omega(iOmega)
              call conjgrad(ham(:,iSpin),&
                  & eiTmp(:nFilled(iSpin)),over, &
                  & ciTmp(:,:nFilled(iSpin),iSpin), &
                  & dci(:,:nFilled(iSpin),jj,iSpin), &
                  & neighborList%iNeighbor, nNeighbor, img2CentCell,iPair, &
                  & iAtomStart, orb%mOrb,tol=tol)
              write(*,*)'Leaving CG'

            end do

            ALLOCATE_(arrayTmp,(size(dci,dim=1),nFilled(iSpin)))
            ALLOCATE_(arrayTmp2,(size(dci,dim=1),size(dci,dim=1)))

            do jj = 1, 2
              arrayTmp = 0.0_dp
              arrayTmp2 = 0.0_dp
              call gemm(arrayTmp,PcT(:,:,iSpin),dci(:,:nFilled(iSpin),jj,iSpin))

              dci(:,:nFilled(iSpin),jj,iSpin) = arrayTmp

              !dci(:,:nFilled(1)) = matmul(PcT,dci(:,:nFilled(1)))

              !dci = ( matmul(dci(:,:nFilled(1)), &
              !    & transpose(ci(:,:nFilled(1),1))) + &
              !    & matmul(ci(:,:nFilled(1),1), &
              !    & transpose(dci(:,:nFilled(1)))) )

              arrayTmp = ci(:,:nFilled(iSpin),iSpin)

              call gemm(arrayTmp2,dci(:,:nFilled(iSpin),jj,iSpin), &
                  & arrayTmp,transB='T')

              arrayTmp = dci(:,:nFilled(iSpin),jj,iSpin)

              call gemm(arrayTmp2,ci(:,:nFilled(iSpin),iSpin),arrayTmp, &
                  & transB='T',beta=1.0_dp)
              arrayTmp2 = 0.25_dp * (arrayTmp2 + transpose(arrayTmp2))
              call packHS(drho(:,iSpin), arrayTmp2,neighborlist%iNeighbor,&
                  & nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
            end do

            DEALLOCATE_(arrayTmp)
            DEALLOCATE_(arrayTmp2)

          end do

          if (nSpin == 1) then
            drho = 2.0_dp * drho
          end if

          call ud2qm(drho)

          dqOut = 0.0_dp
          do iSpin = 1, nSpin
            call mulliken(dqOut(:,:,iSpin), over, drho(:,iSpin), orb, &
                & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
            if (tDFTBU) then
              qBlockOut(:,:,:,iSpin) = 0.0_dp
              call mulliken(qBlockOut(:,:,:,iSpin), over, drho(:,iSpin), &
                  &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
            end if
          end do

          if (tSCC) then
            qOutRed = 0.0_dp
            call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
                & qOutRed(1:nIneqOrb))
            if (tDFTBU) then
              call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
                  & qOutRed )
            end if

            qDiffRed(:) = qOutRed(:) - qInpRed(:)
            sccErrorQ = maxval(abs(qDiffRed))

            write(*,*)iSCCIter,sccErrorQ
            tConverged = (sccErrorQ < sccTol)

            if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
              if (iSCCIter == 1 .and. iOmega == 1) then
                dqIn(:,:,:) = dqOut(:,:,:)
                qInpRed(:) = qOutRed(:)
                if (tDFTBU) then
                  qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
                end if
              else
                call mix(pChrgMixer, qInpRed, qDiffRed)

                call OrbitalEquiv_expand(qInpRed(1:nIneqOrb), iEqOrbitals, &
                    & orb, dqIn)
                if (tDFTBU) then
                  qBlockIn = 0.0_dp
                  call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
                      & qBlockIn, specie0, nUJ, niUJ, iUJ, &
                      & orbEquiv=iEqOrbitals )
                end if
              end if
            end if
          else
            tConverged = .true.
          end if

          if (tConverged) then
            exit lpSCC
          end if

          if (tSpin) then
            chargePerShell(:,:,:) = 0.0_dp
            do iAtom = 1, nAtom
              iSp1 = species(iAtom)
              do iSh1 = 1, orb%nShell(iSp1)
                chargePerShell(iSh1,iAtom,1:nSpin) = &
                    & chargePerShell(iSh1,iAtom,1:nSpin) + &
                    & sum(dqIn(orb%posShell(iSh1,iSp1): &
                    & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
              end do
            end do

          end if

          iSCCIter = iSCCIter +1

        end do lpSCC


        do ii = 1, 3
          write(*,"(E20.12)",advance ='no')-sum(sum(dqOut(:,:,1),dim=1)&
              &*coord0(ii&
              &,:))
          polarisability(ii,iCart) = -sum(sum(dqOut(:,:,1),dim=1)*coord0(ii,:))
        end do
        write(*,*)
        write(1001+icart,"(5E20.12)")omega(iOmega)*Hartree__eV, &
            & polarisability(:,iCart), sqrt(sum((dqOut(:,:,1)*dqOut(:,:,1))))
        !do iAtom = 1, nAtom
        !  write(*,*)dqOut(:,iAtom,1)
        !end do
        !write(*,*)

      end do

    end do

    call destroy(dpotential)
    DEALLOCATE_(eiTmp)
    DEALLOCATE_(ciTmp)
    DEALLOCATE_(dci)
    DEALLOCATE_(dei)
    DEALLOCATE_(hprime)

    DEALLOCATE_(HprimeSqr)
    DEALLOCATE_(SSqr)
    DEALLOCATE_(Pc)
    DEALLOCATE_(PcT)

    if (tWriteTagged) then
      call writeTagged(fdTagged, tag_StatEPol, polarisability)
      close(fdTagged)
    end if

  end subroutine sternheimerDynamic2

  ! sparse but more so
  subroutine sternheimerDynamic3(filling, ei, ci, ham, species, omega, tol, &
      & tSparseBlas)
    real(dp), intent(in)            :: filling(:,:,:)
    real(dp), intent(in)            :: ei(:,:,:)
    real(dp), intent(inout)         :: ci(:,:,:)
    real(dp), intent(in)            :: ham(:,:)
    integer, intent(in)             :: species(:)
    real(dp), intent(in)            :: omega(:)
    real(dp), intent(in)            :: tol
    logical, intent(in)             :: tSparseBlas

    integer  :: ii, jj, kk, iCart, iAtom
    integer  :: nFilled(nSpin)
    integer  :: iSpin, iSpin2, iSh1, iSp1
    real(dp), allocatable :: hprime(:,:)
    type(TPotentials) :: dpotential
    real(dp), allocatable :: HprimeSqr(:,:)
    real(dp), allocatable :: SSqr(:,:)
    real(dp), allocatable :: dci(:,:,:,:)
    real(dp) :: drho(size(over),nSpin)
    real(dp) :: dqIn(orb%mOrb,nAtom,nSpin)
    real(dp) :: dqOut(orb%mOrb,nAtom,nSpin)
    real(dp), allocatable :: arrayTmp(:,:)
    real(dp), allocatable :: ciTmp(:,:,:), eiTmp(:)

    real(dp) :: sccErrorQ, polarisability(3,3)
    integer  :: iSCCIter, iOmega, nOmega
    logical  :: tStopSCC, tConverged, tMetallic

    logical, parameter :: tAnalyse = .false.
    real(dp), allocatable :: dciMag(:,:), dciDecomp(:,:)
    integer :: dcLocal(2)

    nOmega = size(omega)

    if (tWriteTagged) then
      open(fdTagged, file=taggedOut, position="append")
      call writeTagged(fdTagged, tag_StatPol, .true.)
    end if

    if (.not.tSparseBlas) then
      ALLOCATE_(HprimeSqr,(size(ci,dim=1),size(ci,dim=2)))
    else
      ALLOCATE_(HprimeSqr,(0,0))
    end if
    ALLOCATE_(SSqr,(size(ci,dim=1),size(ci,dim=2)))

    do iSpin = 1, nSpin
      nFilled(iSpin) = floor(nEl(iSpin)/real(3-nSpin,dp))
    end do

    call create(dpotential,orb,nAtom,nSpin)

    tMetallic = .false.
    do iSpin = 1, nSpin
      if ( any(filling(:nFilled(ispin),1,iSpin) /= real(3-nSpin,dp)) .or. &
          & any(filling(nFilled(iSpin)+1:,1,iSpin) /= 0.0_dp) ) then
        tMetallic = .true.
        call error("Apparently metallic in static Sternheimer routine, &
            &not implemented yet.")
      end if
    end do

    ALLOCATE_(hprime,(size(over),nSpin))
    ALLOCATE_(dci,(size(ci,dim=1),size(ci,dim=2),2,nSpin))
    ALLOCATE_(arrayTmp,(size(ci,dim=1),size(ci,dim=1)))

    ALLOCATE_(ciTmp,(size(ci,dim=1),maxval(nFilled),nSpin))
    ALLOCATE_(eiTmp,(maxval(nFilled)))
    eiTmp = 0.0_dp

    do iCart = 1, 3

      write(*,*)'Polarisability along ',direction(iCart)

      dci = 1.0_dp
      dqIn = 0.0_dp

      dpotential%extAtom = 0.0_dp
      do iAtom = 1, nAtom
        dpotential%extAtom(iAtom,1) = coord0(iCart,iAtom)
      end do

      dpotential%extBlock = 0.0_dp
      call copy_shiftTo(dpotential%extBlock,dpotential%extAtom, &
          & orb,species)

      do iOmega = 1, nOmega

        dqOut = 0.0_dp


        write(*,*)'At',omega(iOmega),'au',omega(iOmega)*Hartree__eV,'eV'

        !! (Re)Initialize mixer
        if (tSCC) then
          call reset(pChrgMixer,nMixElements) !nIneqOrb)
          if (iOmega == 1) then
            chargePerShell = 0.0_dp
            qBlockIn = 0.0_dp
            qInpRed = 0.0_dp
          end if
        end if
        iSCCIter = 1
        tStopSCC = .false.

        lpSCC: do while (iSCCiter <= nSCCIter)

          dpotential%intAtom = 0.0_dp
          dpotential%intBlock = 0.0_dp
          dpotential%intShell = 0.0_dp
          dpotential%orbitalBlock = 0.0_dp
          dpotential%extAtom = 0.0_dp

          dpotential%orbitalBlock = dpotential%extBlock
          if (tSCC .and. (iSCCiter>1 .or. iOmega > 1)) then
            call updateCharges_SCC(dqIn, orb, species, &
                &neighborList%iNeighbor, img2CentCell)
            call getShiftPerAtom(dpotential%intAtom)
            call getShiftPerL(dpotential%intShell)

            call total_shift(dpotential%intShell,dpotential%intAtom, &
                & dpotential%intShell,orb,species)

            if (tSpin) then
              call addSpinShift(dpotential%intShell,chargePerShell,species,orb&
                  &,W)
            end if
            if (tDFTBU) then
              call shift_DFTBU(dpotential%orbitalBlock,qBlockIn,specie,orb, &
                  & 2, UJ, nUJ, niUJ, iUJ)
            end if
          end if

          call total_shift(dpotential%intBlock,dpotential%intShell, &
              & dpotential%orbitalBlock,orb,species)

          hprime = 0.0_dp

          call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
              & species,orb,iPair,nAtom,img2CentCell,dpotential%intBlock)

          if (nSpin == 2) then
            hprime(:,:) = 2.0_dp *  hprime(:,:)
          end if
          call qm2ud(hprime)

          !dei = 0.0_dp
          drho = 0.0_dp

          do iSpin = 1, nSpin

            if (.not.tSparseBlas) then ! H' matrix as dense
              HprimeSqr = 0.0_dp
              call unpackHS(HprimeSqr,hprime(:,iSpin),neighborList%iNeighbor, &
                  & nNeighbor, iAtomStart,iPair,img2CentCell)
              do ii = 1, size(ci,dim=2)
                HprimeSqr(ii,ii+1:) = HprimeSqr(ii+1:,ii)
              end do
            end if

            if (tStoreEigvecs) then
              iSpin2 = 1
              call get(storeEigvecsReal(iSpin), ci(:,:,iSpin2))
            else
              iSpin2 = iSpin
            end if

            !dei(1:nFilled(iSpin),iSpin) = sum(ci(:,1:nFilled(iSpin),iSpin2) &
            !    & *matmul(HprimeSqr,ci(:,1:nFilled(iSpin),iSpin2)),dim=1)

            if (tSparseBlas) then
              call dftbSYMM(ciTmp(:,:nFilled(iSpin),iSpin),hprime(:,iSpin), &
                  & ci(:,:nFilled(iSpin),iSpin2),neighborList%iNeighbor, &
                  & nNeighbor, img2CentCell, iPair, iAtomStart, orb%mOrb, &
                  & alpha=-1.0_dp)
            else
              ciTmp(:,:nFilled(iSpin),iSpin) = -matmul(HprimeSqr, &
                  & ci(:,:nFilled(iSpin),iSpin2))
              !call gemm(ciTmp(:,:nFilled(iSpin),iSpin),HprimeSqr, &
              !    & ci(:,:nFilled(iSpin),iSpin2),alpha=-1.0_dp)
            end if

            do jj = 1, 2
              kk = 2*jj-3

              write(*,*)'entering sparse CG routine'
              eiTmp(:nFilled(iSpin)) = &
                  & ei(:nFilled(iSpin),1,iSpin)+kk*omega(iOmega)
              call conjgrad(ham(:,iSpin),&
                  & eiTmp(:nFilled(iSpin)),over, &
                  & ciTmp(:,:nFilled(iSpin),iSpin), &
                  & dci(:,:nFilled(iSpin),jj,iSpin), &
                  & neighborList%iNeighbor, nNeighbor, img2CentCell,iPair, &
                  & iAtomStart, orb%mOrb,tol=tol)
              write(*,*)'Leaving CG routine'

            end do

            do jj = 1, 2
              arrayTmp = 0.0_dp

              arrayTmp = matmul(dci(:,:nFilled(iSpin),jj,iSpin), &
                  & transpose(ci(:,:nFilled(iSpin),iSpin)))
              !call gemm(arrayTmp,dci(:,:nFilled(iSpin),jj,iSpin), &
              !    & ci(:,:nFilled(iSpin),iSpin),transB='T')
              arrayTmp = 0.5_dp * (arrayTmp + transpose(arrayTmp))
              call packHS(drho(:,iSpin), arrayTmp,neighborlist%iNeighbor,&
                  & nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)

            end do

          end do

          if (nSpin == 1) then
            drho = 2.0_dp * drho
          end if

          call ud2qm(drho)

          dqOut = 0.0_dp
          do iSpin = 1, nSpin
            call mulliken(dqOut(:,:,iSpin), over, drho(:,iSpin), orb, &
                & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
            if (tDFTBU) then
              qBlockOut(:,:,:,iSpin) = 0.0_dp
              call mulliken(qBlockOut(:,:,:,iSpin), over, drho(:,iSpin), &
                  &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
            end if
          end do

          if (tSCC) then
            qOutRed = 0.0_dp
            call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
                & qOutRed(1:nIneqOrb))
            if (tDFTBU) then
              call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
                  & qOutRed )
            end if

            qDiffRed(:) = qOutRed(:) - qInpRed(:)
            sccErrorQ = maxval(abs(qDiffRed))

            write(*,*)iSCCIter,sccErrorQ
            tConverged = (sccErrorQ < sccTol)

            if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
              if (iSCCIter == 1 .and. iOmega == 1) then
                dqIn(:,:,:) = dqOut(:,:,:)
                qInpRed(:) = qOutRed(:)
                if (tDFTBU) then
                  qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
                end if
              else
                call mix(pChrgMixer, qInpRed, qDiffRed)

                call OrbitalEquiv_expand(qInpRed(1:nIneqOrb), iEqOrbitals, &
                    & orb, dqIn)
                if (tDFTBU) then
                  qBlockIn = 0.0_dp
                  call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
                      & qBlockIn, specie0, nUJ, niUJ, iUJ, &
                      & orbEquiv=iEqOrbitals )
                end if
              end if
            end if
          else
            tConverged = .true.
          end if

          if (tConverged) then
            exit lpSCC
          end if

          if (tSpin) then
            chargePerShell(:,:,:) = 0.0_dp
            do iAtom = 1, nAtom
              iSp1 = species(iAtom)
              do iSh1 = 1, orb%nShell(iSp1)
                chargePerShell(iSh1,iAtom,1:nSpin) = &
                    & chargePerShell(iSh1,iAtom,1:nSpin) + &
                    & sum(dqIn(orb%posShell(iSh1,iSp1): &
                    & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
              end do
            end do

          end if

          iSCCIter = iSCCIter +1

        end do lpSCC

        if (tAnalyse) then
          ALLOCATE_(dciMag,(maxval(nFilled),2))
          ALLOCATE_(dciDecomp,(size(ci,dim=2),2))
          write(*,*)'Magnitude of KS state changes'
          do iSpin = 1, nSpin
            dciMag = 0.0_dp
            do jj = 1, 2
              call dftbSYMM( ciTmp(:,:nFilled(iSpin),iSpin), over, &
                  & dci(:,:nFilled(iSpin),jj,iSpin),neighborList%iNeighbor, &
                  & nNeighbor, img2CentCell, iPair, iAtomStart, orb%mOrb )
              do ii = 1,nFilled(iSpin)
                dciMag(ii,jj)=sqrt( dot_product(ciTmp(:,ii,iSpin), &
                    & dci(:,ii,jj,iSpin)) )
              end do
            end do
            dciDecomp = 0.0_dp
            do jj = 1, 2
              dcLocal(jj:jj) = maxloc(dciMag(:,jj))
              call dftbSYMv(ciTmp(:,1,1),over,dci(:,dcLocal(jj),jj,iSpin), &
                  & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, &
                  & iAtomStart, orb%mOrb)
              do ii = nFilled(iSpin)+1, size(ci,dim=2)
                dciDecomp(ii,jj) = dot_product( ciTmp(:,1,1), ci(:,ii,iSpin) )
              end do
            end do
            write(*,*)'Largest change in ground state levels',dcLocal
            do jj = 1, 2
              dcLocal(jj:jj) = maxloc(abs(dciDecomp(:,jj)))
            end do
            write(*,*)'Transition to',dcLocal
            do ii = nFilled(iSpin)+1, size(ci,dim=2)
              write(*,*)ii,dciDecomp(ii,:)
            end do
          end do
          DEALLOCATE_(dciMag)
          DEALLOCATE_(dciDecomp)
          write(*,*)
        end if

        write(*,"(1X,A,1X)",advance ='no')'mu'
        do ii = 1, 3
          write(*,"(E20.12)",advance ='no') &
              & -sum(sum(dqOut(:,:,1),dim=1)*coord0(ii,:))
          polarisability(ii,iCart) = -sum(sum(dqOut(:,:,1),dim=1)*coord0(ii,:))
        end do
        write(*,*)
        write(1001+icart,"(5E20.12)")omega(iOmega)*Hartree__eV, &
            & polarisability(:,iCart), sqrt(sum((dqOut(:,:,1)*dqOut(:,:,1))))

      end do

    end do

    call destroy(dpotential)
    DEALLOCATE_(eiTmp)
    DEALLOCATE_(ciTmp)
    DEALLOCATE_(dci)
    DEALLOCATE_(hprime)

    DEALLOCATE_(HprimeSqr)
    DEALLOCATE_(SSqr)
    DEALLOCATE_(arrayTmp)

    if (tWriteTagged) then
      call writeTagged(fdTagged, tag_StatEPol, polarisability)
      close(fdTagged)
    end if

  end subroutine sternheimerDynamic3
  
end module dftbp_derivs_sternheimer
