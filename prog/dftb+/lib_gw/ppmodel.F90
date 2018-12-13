!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to calculate the correlation part of the GW self energy using a plasmon pole
!> model
module PPModel
  use accuracy, only : dp
  use eigensolver
  use TransDens
  implicit none

  private
  public plasmonPole, selfCor

contains
  !> Implements a global diagonalizable Plasmon-Pole-Model for [W-v].
  !> This quantity is symmetric in contrast to [eps].
  !> The two test frequencies of the model are pure imaginary.
  !> Ansatz: [W-v](w) = phi*L(w)*phi^T; with the eigenvectors phi of [W-v](w=0)
  !> and L(w) = zq*wq^2/(w^2-wq^2)
  subroutine plasmonPole(nAng, wMinV, omegaPPi, eVecPP, ampPP, frqPP, mode)

    !> Size of GW basis
    integer,          intent(in)    :: nAng

    !> mode=1 : [W-v] is diagonalized, zq determined; mode=2 : wq is determined
    integer,          intent(in)    :: mode

    !> [W-v](w) (overwritten on output)
    real(dp),         intent(inout) :: wMinV(:,:)

    !> Imaginary part of test frequency
    real(dp),         intent(in)    :: omegaPPi

    !> phi
    real(dp),         intent(out)   :: eVecPP(:,:)

    !> zq
    real(dp),         intent(out)   :: ampPP(:)

    !> wq
    complex(dp),      intent(out)   :: frqPP(:)

    integer                          :: iAng
    real(dp)                         :: EigenVal(nAng)
    real(dp)                         :: rTmp

    @:ASSERT(all(shape(wMinV) == (/ nAng, nAng /)))
    @:ASSERT(all(shape(eVecPP) == (/ nAng, nAng /)))
    @:ASSERT(size(frqPP) == nAng)
    @:ASSERT(size(ampPP) == nAng)

    if(mode == 1) then
      @:ASSERT(abs(omegaPPi) < epsilon(1.0_dp))
      call heev(wMinV, EigenVal, 'L', 'V')
      eVecPP(:,:) = wMinV(:,:)
      ampPP(:) = -EigenVal(:)
    else
      call heev(wMinV, EigenVal, 'L', 'N')
      do iAng = 1, nAng
        rTmp =  -EigenVal(iAng)*omegaPPi**2/(ampPP(iAng)+EigenVal(iAng))
        if(rTmp < 0.0_dp) then
          frqPP(iAng) = (0.0_dp,1.0_dp)*sqrt(-rTmp)
        else
          frqPP(iAng) = (1.0_dp,0.0_dp)*sqrt(rTmp)
        endif
      enddo
    endif

  end subroutine plasmonPole

  !> Computes the correlation part of the self energy <psi_i|Sigma_c(E)|psi_i> based on a
  !> Plasmon-Pole-Model (see plasmonPole routine)
  subroutine selfCor(nOrb, nAng, nAtom, mAngAtom, eigenVec, filling, eVecPP, ampPP, frqPP,&
      & eigenVal, eVar, eSelfCor)

    !> Size of DFTB basis
    integer,          intent(in)    :: nOrb

    !> Size of GW basis
    integer,          intent(in)    :: nAng

    !> Number of Atoms
    integer,          intent(in)    :: nAtom

    !> Max. angular momentum per atom
    integer,          intent(in)    :: mAngAtom(:)

    !> Eigenvectors of DFTB Hamiltonian
    real(dp),         intent(in)    :: eigenVec(:,:)

    !> Electron occupation of orbital (spin restricted)
    real(dp),         intent(in)    :: filling(:)

    !> Eigenvectors of [W-v]
    real(dp),         intent(in)    :: eVecPP(:,:)

    !> Amplitudes of the PPM
    real(dp),         intent(in)    :: ampPP(:)

    !> Eigenvalues of DFTB Hamiltonian
    real(dp),         intent(in)    :: eigenVal(:)

    !> Energy variable (E from above)
    real(dp),         intent(in)    :: eVar(:)

    !> Frequencies of the PPM
    complex(dp),      intent(in)    :: frqPP(:)

    !> <psi_i|Sigma_c(E)|psi_i>
    complex(dp),      intent(out)   :: eSelfCor(:)

    integer                          :: iOrb, kOrb, iHOMO, iAng
    real(dp), allocatable            :: qTrans(:)
    real(dp)                         :: qTimPhi

    @:ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    @:ASSERT(size(filling) == nOrb)
    @:ASSERT(all(shape(eVecPP) == (/ nAng, nAng /)))
    @:ASSERT(size(frqPP) == nAng)
    @:ASSERT(size(ampPP) == nAng)
    @:ASSERT(size(eigenVal) == nOrb)
    @:ASSERT(size(eVar) == nOrb)
    @:ASSERT(size(eSelfCor) == nOrb)

    allocate(qTrans(nAng))

    ! Find highest occupied orbital
    do iOrb = 1, nOrb
      if(abs(filling(iOrb) < epsilon(1.0_dp))) exit
    enddo
    iHOMO = iOrb-1

    eSelfCor(:) = cmplx(0,0,dp)
    do iOrb = 1, nOrb
      do kOrb = 1, iHOMO
        call transQ(nAtom, mAngAtom, iOrb, kOrb, eigenVec, qTrans)
        do iAng = 1, nAng
          if(abs(ampPP(iAng)) < epsilon(1.0_dp)) exit
          qTimPhi = dot_product(qTrans, eVecPP(:,iAng))
          eSelfCor(iOrb) = eSelfCor(iOrb) + 0.5_dp * filling(kOrb)*&
              & qTimPhi * qTimPhi * (0.5_dp * ampPP(iAng) * frqPP(iAng))/&
              &(eVar(iOrb) - eigenVal(kOrb) + frqPP(iAng))
        enddo
      enddo
      do kOrb = iHOMO+1, nOrb
        call transQ(nAtom, mAngAtom, iOrb, kOrb, eigenVec, qTrans)
        do iAng = 1, nAng
          if(abs(ampPP(iAng)) < epsilon(1.0_dp)) exit
          qTimPhi = dot_product(qTrans, eVecPP(:,iAng))
          eSelfCor(iOrb) = eSelfCor(iOrb) + &
              & qTimPhi * qTimPhi * (0.5_dp * ampPP(iAng) * frqPP(iAng))/&
              &(eVar(iOrb) - eigenVal(kOrb) - frqPP(iAng))
        enddo
      enddo
    enddo

    deallocate(qTrans)

  end subroutine selfCor

end module PPModel
