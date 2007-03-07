!!* Contains subroutines to calculate the correlation part of the GW
!!* self energy using a plasmon pole model
module PPModel
# include "assert.h"
# include "allocate.h"
  use accuracy, only : dp
  use eigensolver
  use TransDens
  implicit none

  private
  public plasmonPole, selfCor  
  
contains
  !!* Implements a global diagonalizable Plasmon-Pole-Model for [W-v].
  !!* This quantity is symmetric in contrast to [eps]. The two test
  !!* frequencies of the model are pure imaginary.
  !!* Ansatz: [W-v](w) = phi*L(w)*phi^T; with the eigenvectors phi of [W-v](w=0)
  !!*         and L(w) = zq*wq^2/(w^2-wq^2)
  !!* @param nAng                   Size of GW basis
  !!* @param wMinV(nAng,nAng)       [W-v](w) (overwritten on output)
  !!* @param omegaPPi               Imaginary part of test frequency
  !!* @param eVecPP(nAng,nAng)      phi
  !!* @param ampPP(nAng)            zq
  !!* @param frqPP(nAng)            wq
  !!* @param mode                   =1 : [W-v] is diagonalized, zq determined
  !!*                               =2 : wq is determined
  subroutine plasmonPole(nAng, wMinV, omegaPPi, eVecPP, ampPP, frqPP, mode)
    integer,          intent(in)    :: nAng
    integer,          intent(in)    :: mode
    real(dp),         intent(inout) :: wMinV(:,:)
    real(dp),         intent(in)    :: omegaPPi
    real(dp),         intent(out)   :: eVecPP(:,:)
    real(dp),         intent(out)   :: ampPP(:)
    complex(dp),      intent(out)   :: frqPP(:)

    integer                          :: iAng
    real(dp)                         :: EigenVal(nAng)
    real(dp)                         :: rTmp

    ASSERT(all(shape(wMinV) == (/ nAng, nAng /)))
    ASSERT(all(shape(eVecPP) == (/ nAng, nAng /)))
    ASSERT(size(frqPP) == nAng)
    ASSERT(size(ampPP) == nAng)

    if(mode == 1) then
      ASSERT(abs(omegaPPi) < epsilon(1.0_dp))
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
  
  !!* Computes the correlation part of the self energy <psi_i|Sigma_c(E)|psi_i>
  !!* based on a Plasmon-Pole-Model (see plasmonPole)
  !!* @param nOrb                 Size of DFTB basis
  !!* @param nAng                 Size of GW basis
  !!* @param nAtom                Number of Atoms
  !!* @param mAngAtom(nAtom)      Max. angular momentum per atom
  !!* @param eigenVec(nOrb,nOrb)  Eigenvectors of DFTB Hamiltonian
  !!* @param filling(nOrb)        Electron occupation of orbital (spin resticted)  
  !!* @param eVecPP(nAng,nAng)    Eigenvectors of [W-v] 
  !!* @param ampPP(nAng)          Amplitudes of the PPM
  !!* @param frqPP(nAng)          Frequencies of the PPM
  !!* @param eigenVal(nOrb)       Eigenvalues of DFTB Hamiltonian  
  !!* @param eVar(nOrb)           Energy variable (E from above)
  !!* @param eSelfCor(nOrb)       <psi_i|Sigma_c(E)|psi_i>                      
  subroutine selfCor(nOrb, nAng, nAtom, mAngAtom, eigenVec, filling, &
      &   eVecPP, ampPP, frqPP, eigenVal, eVar, eSelfCor)
    integer,          intent(in)    :: nOrb
    integer,          intent(in)    :: nAng
    integer,          intent(in)    :: nAtom
    integer,          intent(in)    :: mAngAtom(:)
    real(dp),         intent(in)    :: eigenVec(:,:)
    real(dp),         intent(in)    :: filling(:)
    real(dp),         intent(in)    :: eVecPP(:,:)
    real(dp),         intent(in)    :: ampPP(:)
    real(dp),         intent(in)    :: eigenVal(:)
    real(dp),         intent(in)    :: eVar(:)
    complex(dp),      intent(in)    :: frqPP(:)
    complex(dp),      intent(out)   :: eSelfCor(:)

    integer                          :: iOrb, kOrb, iHOMO, iAng
    real(dp), allocatable            :: qTrans(:)
    real(dp)                         :: qTimPhi

    ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    ASSERT(size(filling) == nOrb)
    ASSERT(all(shape(eVecPP) == (/ nAng, nAng /)))
    ASSERT(size(frqPP) == nAng)
    ASSERT(size(ampPP) == nAng)
    ASSERT(size(eigenVal) == nOrb)
    ASSERT(size(eVar) == nOrb)
    ASSERT(size(eSelfCor) == nOrb)
    ALLOCATE_(qTrans, (nAng))
    
    !! Find highest occupied orbital
    do iOrb = 1, nOrb
      if(abs(filling(iOrb) < epsilon(1.0_dp))) exit
    enddo
    iHOMO = iOrb-1

    eSelfCor(:) = (0.0_dp, 0.0_dp)
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
    
    DEALLOCATE_(qTrans)
  end subroutine selfCor
end module PPModel
