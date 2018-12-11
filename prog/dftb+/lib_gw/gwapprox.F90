!!* Implements the GW approximation for DFTB 
module GWApprox
#include "allocate.h"  
#include "assert.h"
  use accuracy, only : dp, deltaSelfEnergy
  use blasroutines
  use TransDens
  use GammaMat
  use RandomPhase
  use PPModel
  use gwXandXC
  implicit none
  private

  public :: gwDriver
  
contains
  !!* Implements the GW approximation for DFTB based on the article
  !!* ==> Niehaus et. al., PRA 71, 022508 (2005)
  !!* @param nAtom                Number of Atoms
  !!* @param nType                Number of different species
  !!* @param mAngAtom(nAtom)      Max. angular momentum per atom
  !!* @param mAngSpecie(nType)    Max. angular momentum per species
  !!* @param iAtomStart(nAtom)    Indexing array for square Hamiltonian
  !!* @param specie(nAtom)        Species list of atomic species for all atom
  !!* @param coord(3,nAtom)       Atomic coordinates
  !!* @param eigenVec(nOrb,nOrb)  Eigenvectors of DFTB Hamiltonian
  !!* @param eigenVal(nOrb,nOrb)  Eigenvalues of DFTB Hamiltonian
  !!* @param SSqrReal(nOrb,nOrb)  Overlap matrix (square format)
  !!* @param XCSqr(nOrb,nOrb)     Matrix of Exchange-Correlation potential
  !!* @param filling(nOrb)        Electron occupation of orbital (spin restr.)
  !!* @param eeU(l,nType)         Electron repulsion integrals (compressed WF)
  !!* @param xTab(nType,m,m)      Atomic exchange integrals
  !!* @param eTab(nType,m,m)      Electron repulsion integrals (uncompressed WF)
  !!* @param hubbU(l,nType)       Hubbard values (ee + xc)
  !!* @param qOutput(m,nAtom)     m-resolved populations
  !!* @param referenceN0(l,nType) l-resolved atomic populations 
  !!* @param eHFX(nOrb)           Hartree-Fock exchange energy per orb.
  !!* @param eXC(nOrb)            Exchange-Correlation energy per orb.
  !!* @param eQP(nOrb)            Quasi particle energies
  subroutine gwDriver(nAtom, nType, mAngAtom, mAngSpecie, iAtomStart,&
      & specie, coord, eigenVec, eigenVal, SSqrReal, XCSqr,          & 
      & filling , eeU, xTab, eTab, hubbU, qOutput, referenceN0,      &
      & eHFX, eXC, eQP)
    
    integer,  intent(in)   :: nAtom
    integer,  intent(in)   :: nType
    integer,  intent(in)   :: mAngAtom(:)
    integer,  intent(in)   :: mAngSpecie(:)
    integer,  intent(in)   :: iAtomStart(:)
    integer,  intent(in)   :: specie(:)
    real(dp), intent(in)   :: coord(:,:)
    real(dp), intent(in)   :: eigenVec(:,:)
    real(dp), intent(in)   :: eigenVal(:)
    real(dp), intent(in)   :: SSqrReal(:,:)
    real(dp), intent(in)   :: XCSqr(:,:)
    real(dp), intent(in)   :: filling(:)
    real(dp), intent(in)   :: eeU(:,:)
    real(dp), intent(in)   :: xTab(:,:,:)
    real(dp), intent(in)   :: eTab(:,:,:)
    real(dp), intent(in)   :: hubbU(:,:)
    real(dp), intent(in)   :: qOutput(:,:)
    real(dp), intent(in)   :: referenceN0(:,:)
    real(dp), intent(out)  :: eHFX(:)
    real(dp), intent(out)  :: eXC(:)
    real(dp), intent(out)  :: eQP(:)

    integer                  :: nAng, maxM
    integer                  :: iOrb,jOrb,nOrb
    integer                  :: iPPM
    real(dp), allocatable    :: qTrans(:)
    real(dp), allocatable    :: delEigenVal(:,:)
    real(dp), allocatable    :: eeGamma(:,:)
    real(dp), allocatable    :: epsInv(:,:)
    real(dp), allocatable    :: wMinV(:,:)
    real(dp), allocatable    :: eVecPP(:,:)
    real(dp), allocatable    :: ampPP(:)
    real(dp), allocatable    :: eVar(:)
    real(dp), parameter      :: omegaPPi(2)=(/0.0_dp,0.5_dp/)
    complex(dp), allocatable :: frqPP(:)
    complex(dp), allocatable :: eSelfCor(:)
    complex(dp), allocatable :: zRenorm(:)
    complex(dp), allocatable :: delQP(:)
    
    
    nAng = sum(mAngAtom) + nAtom
    nOrb = dot_product(mAngAtom+1,mAngAtom+1)
    ASSERT(size(mAngAtom) == nAtom)
    ASSERT(size(mAngSpecie) == nType)
    ASSERT(size(iAtomStart) == nAtom)
    ASSERT(size(specie) == nAtom)
    ASSERT(all(shape(coord) == (/ 3, nAtom /)))
    ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    ASSERT(size(eigenVal) == nOrb)
    ASSERT(all(shape(SSqrReal) == shape(eigenVec)))
    ASSERT(all(shape(XCSqr) == shape(eigenVec)))
    ASSERT(size(filling) == nOrb)
    ASSERT(all(shape(eeU) == (/ maxval(mAngAtom+1), nType /)))
    ASSERT(all(shape(hubbU) == (/ maxval(mAngAtom+1), nType /)))
    maxM = (maxval(mAngAtom)+1)**2
    ASSERT(all(shape(eTab) == (/ nType, maxM, maxM /)))
    ASSERT(all(shape(xTab) == (/ nType, maxM, maxM /)))
    ASSERT(all(shape(qOutput) == (/ maxval(mAngAtom+1)**2, nAtom /)))
    ASSERT(all(shape(referenceN0) == (/ maxval(mAngAtom+1), nType /)))
    ASSERT(size(eHFX) == nOrb)
    ASSERT(size(eXC) == nOrb)
    ASSERT(size(eQP) == nOrb)
    
    ALLOCATE_(qTrans, (nAng))
    call transQ_init(SSqrReal,eigenVec)
    ALLOCATE_(delEigenVal, (nOrb,nOrb))
    ALLOCATE_(eeGamma, (nAng,nAng))
    ALLOCATE_(epsInv, (nAng,nAng))
    ALLOCATE_(wMinV, (nAng,nAng))
    ALLOCATE_(eVecPP, (nAng,nAng))
    ALLOCATE_(ampPP, (nAng))
    ALLOCATE_(frqPP, (nAng))    
    ALLOCATE_(eVar, (nOrb))
    ALLOCATE_(eSelfCor, (nOrb))
    ALLOCATE_(zRenorm, (nOrb))
    ALLOCATE_(delQP, (nOrb))

    print *,'       ==> Start of GW calculation'
    do jOrb = 1,nOrb
      do iOrb = 1,jOrb-1
        delEigenVal(iOrb,jOrb) = eigenVal(jOrb) - eigenVal(iOrb)
      enddo
    enddo
    !! Get matrix of the coulomb interaction [v] (equal to gamma)
    !! eeU Hubbard parameter contain no XC contribution
    call gammaM(nAtom, mAngAtom, specie, coord, eeU, eeGamma)

    !! Plasmon-Pole-Model with two test frequencies omegaPPi
    do iPPM = 1,2
      !! ([eps(w)]S)^-1 = (1 - [v]<P(w)>)^-1
      call buildEpsInv(nAtom, nOrb, nAng, mAngAtom, filling, &
          & eigenVec, delEigenVal, eeGamma, omegaPPi(iPPM), epsInv)
      !! wMinV =  epsInv * [v] - [v]
      wMinV(:,:) = eeGamma(:,:)
      call symm(wMinV, 'R', eeGamma, epsInv, 'L', 1.0_dp, -1.0_dp)
      call plasmonPole(nAng, wMinV, omegaPPi(iPPM), eVecPP, &
          & ampPP, frqPP, iPPM)
    enddo
    !! PPM parameters plasmon frequencies (frqPP) and amplitudes (ampPP)
    !! are now known. eVecPP are eigenvectors of wMinV

    !! Get Hartree-Fock exchange (eHFX) and v_xc (eXC) energies
    !! Implements Eq. 17 of the PRA paper
    call hfExchange(nOrb, nAng, nAtom,  mAngAtom, iAtomStart, &
        & specie, eeGamma, eigenVec, filling, xTab, eHFX)
    !! Implements Eq. 22b 
    call dftbXC(nOrb, nAng, nAtom, mAngAtom, mAngSpecie, specie, &
        & coord, XCSqr, eigenVec, eTab, hubbU, qOutput, referenceN0, eXC)

    !! Consider energy dependence of self energy through renormalization
    !! Self energy is numerically differentiated (Eq. 3) 
    eVar(:) = eigenVal(:) + deltaSelfEnergy
    !! Computation of the correlation part of the self energy (Eq. 21)
    call selfCor(nOrb, nAng, nAtom, mAngAtom, eigenVec, filling, &
        &   eVecPP, ampPP, frqPP, eigenVal, eVar, eSelfCor)
    zRenorm(:) = eSelfCor(:)
    eVar(:) = eigenVal(:)
    call selfCor(nOrb, nAng, nAtom, mAngAtom, eigenVec, filling, &
        &   eVecPP, ampPP, frqPP, eigenVal, eVar, eSelfCor)
    zRenorm(:) = (zRenorm(:) -  eSelfCor(:))/deltaSelfEnergy
    zRenorm(:) = 1.0_dp/(1.0_dp - zRenorm(:))
    
    !! Final quasi particle energies (Eq. 2)
    delQP(:) = zRenorm(:)*(eSelfCor(:) + eHFX(:) - eXC(:))
    eQP(:) = eigenVal(:) + real(delQP(:))

    DEALLOCATE_(delQP)
    DEALLOCATE_(zRenorm)
    DEALLOCATE_(eSelfCor)
    DEALLOCATE_(eVar)
    DEALLOCATE_(frqPP)    
    DEALLOCATE_(eVecPP)    
    DEALLOCATE_(ampPP)    
    DEALLOCATE_(epsInv)
    DEALLOCATE_(wMinV)
    DEALLOCATE_(eeGamma)
    DEALLOCATE_(delEigenVal)
    call transQ_destr()
    DEALLOCATE_(qTrans)
    print *,'       ==> GW calculation finished'
  end subroutine gwDriver

end module GWApprox
