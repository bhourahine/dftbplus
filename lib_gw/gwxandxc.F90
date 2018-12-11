!!* Contains subroutines to calculate the Hartee-Fock exchange part of the GW
!!* self energy and the DFT Exchange and correlation contribution
module gwXandXC
# include "assert.h"
# include "allocate.h"
  use accuracy, only : dp
  use blasroutines
  use TransDens
  use GammaMat
  implicit none

  private
  public hfExchange, dftbXC
  
contains

  !!* Calculates the HF exchange matrix elements <psi_i|Sigma_x|psi_i>
  !!* in the basis of DFTB Kohn-Sham Orbitals. The exchange is calculated in the
  !!* gamma-approximation with Hubbard parameters derived from Coulomb
  !!* interaction only. Onsite exchange elements are accounted for additionaly.
  !!* @param nOrb                 Size of DFTB basis
  !!* @param nAng                 Size of GW basis
  !!* @param nAtom                Number of Atoms
  !!* @param mAngAtom(nAtom)      Max. angular momentum per atom 
  !!* @param iAtomStart(nAtom)    Indexing array for square Hamiltonian
  !!* @param specie(nAtom)        Species list of atomic species for all atom
  !!* @param eeGamma(nAng,nAng)   Gamma matrix (pure Coulomb)
  !!* @param eigenVec(nOrb,nOrb)  Eigenvectors of Hamiltonian
  !!* @param filling(nOrb)        Electron occupation of orbital (spin resticted)
  !!* @param xTab(nType,:,:)      Atomic exchange integrals
  !!* @param eHFX(nOrb)           <psi_i|Sigma_x|psi_i>
  subroutine hfExchange(nOrb, nAng, nAtom, mAngAtom, iAtomStart, specie, eeGamma,&
                      &  eigenVec, filling, xTab, eHFX)

    integer,          intent(in)    :: nOrb
    integer,          intent(in)    :: nAng
    integer,          intent(in)    :: nAtom
    integer,          intent(in)    :: mAngAtom(:)
    integer,          intent(in)    :: iAtomStart(:)
    integer,          intent(in)    :: specie(:)
    real(dp),         intent(in)    :: eeGamma(:,:)
    real(dp),         intent(in)    :: eigenVec(:,:)
    real(dp),         intent(in)    :: filling(:)
    real(dp),         intent(in)    :: xtab(:,:,:)
    real(dp),         intent(out)   :: eHFX(:)

    integer                          :: nType, maxM 
    integer                          :: iOrb, jOrb, iHOMO, iAtom
    integer                          :: iM1, iM2, ind1, ind2
    real(dp), allocatable            :: qTrans(:)
    real(dp), allocatable            :: vTmp(:)
    real(dp), allocatable            :: mTmp(:,:)
    real(dp)                         :: alpha
    
    ASSERT(size(specie) == nAtom)
    ASSERT(all(shape(eeGamma) == (/ nAng, nAng /)))
    ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    ASSERT(size(filling) == nOrb)
    nType = maxval(specie)
    maxM = (maxval(mAngAtom)+1)**2
    ASSERT(all(shape(xTab) == (/ nType, maxM, maxM /)))
    ASSERT(size(eHFX) == nOrb)
    ALLOCATE_(qTrans, (nAng))
    ALLOCATE_(vTmp, (nAng))
    ALLOCATE_(mTmp, (nOrb, nOrb))
    
    !! Find highest occupied orbital
    do iOrb = 1, nOrb
      if(abs(filling(iOrb) < epsilon(1.0_dp))) exit
    enddo
    iHOMO = iOrb-1

    !! E_X = -sum_ij q_ij gamma q_ij (j\in occ)
    vTmp(:) = 0.0_dp
    eHFX(:) = 0.0_dp
    do iOrb = 1, nOrb
      do jOrb = 1, iHOMO
        call transQ(nAtom, mAngAtom, iOrb, jOrb, eigenVec, qTrans)
        alpha =  -0.5_dp*filling(jOrb)
        call symv(vTmp, eeGamma, qTrans, 'L', alpha, 0.0_dp)
        eHFX(iOrb) = eHFX(iOrb) + dot_product(qTrans, vTmp)
      enddo
    enddo

    !! mtmp = sum_j c(m,j)*c(n,j) (j\in occ) 
    mTmp(:,:) = 0.0_dp
    do jOrb = 1, iHOMO
       call ger(mTmp, 0.5_dp*filling(jOrb), eigenVec(:,jOrb), eigenVec(:,jOrb))
    enddo

    !! Add onsite exchange integrals
    do iOrb = 1, nOrb
      do iAtom = 1, nAtom
        do iM1 = 1, (mAngAtom(iAtom)+1)**2
          ind1 = iAtomStart(iAtom) + iM1 - 1
          do iM2 = iM1+1, (mAngAtom(iAtom)+1)**2
            ind2 = iAtomStart(iAtom) + iM2 - 1
            eHFX(iOrb) = eHFX(iOrb) - xTab(specie(iAtom), iM1, iM2)*&
                &( mTmp(ind2,ind2)*(eigenVec(ind1,iOrb)**2)&
                &+ mTmp(ind1,ind1)*(eigenVec(ind2,iOrb)**2)&
                &+ 2*mTmp(ind1,ind2)*eigenVec(ind2,iOrb)*eigenVec(ind1,iOrb))
          enddo
        enddo
      enddo
    enddo
    
    DEALLOCATE_(mTmp)    
    DEALLOCATE_(vTmp)
    DEALLOCATE_(qTrans)
  end subroutine hfExchange
  
  !!* Calculates the exchange correlation matrix elements <psi_i|v_xc|psi_i> 
  !!* in the basis of DFTB Kohn-Sham Orbitals. They are calculated from xc-only
  !!* SK tables. To account for the density dependence of v_xc, a SCC
  !!* correction is applied.
  !!* @param nOrb                 Size of DFTB basis
  !!* @param nAng                 Size of GW basis
  !!* @param nAtom                Number of Atoms
  !!* @param mAngAtom(nAtom)      Max. angular momentum per atom
  !!* @param mAngSpecie(nType)    Max. angular momentum per species
  !!* @param specie(nAtom)        Species list of atomic species for all atom  
  !!* @param coord(3,nAtom)       Atomic coordinates  
  !!* @param XCSqr(nOrb,nOrb)     Matrix of Exchange-Correlation potential 
  !!* @param eigenVec(nOrb,nOrb)  Eigenvectors of Hamiltonian
  !!* @param eTab(nType,:,:)      Atomic repulsion integrals (uncompressed)
  !!* @param hubbU(l,nType)       Hubbard values (ee + xc)
  !!* @param qOutput(m,nAtom)     m-resolved populations
  !!* @param referenceN0(l,nType) l-resolved atomic populations   
  !!* @param eXC(nOrb)            <psi_i|v_xc|psi_i> 
  subroutine dftbXC(nOrb, nAng, nAtom, mAngAtom, mAngSpecie, specie, coord,  &
      &             XCSqr, eigenVec, eTab, hubbU, qOutput, referenceN0, eXC)

    integer,          intent(in)    :: nOrb
    integer,          intent(in)    :: nAng
    integer,          intent(in)    :: nAtom
    integer,          intent(in)    :: mAngAtom(:)
    integer,          intent(in)    :: mAngSpecie(:)
    integer,          intent(in)    :: specie(:)
    real(dp),         intent(in)    :: coord(:,:)
    real(dp),         intent(in)    :: XCSqr(:,:)
    real(dp),         intent(in)    :: eigenVec(:,:)
    real(dp),         intent(in)    :: etab(:,:,:)
    real(dp),         intent(in)    :: hubbU(:,:)
    real(dp),         intent(in)    :: qOutput(:,:)
    real(dp),         intent(in)    :: referenceN0(:,:)
    real(dp),         intent(out)   :: eXC(:)

    integer                          :: nType, maxM
    integer                          :: iOrb, iType, iM1, iM2, mmAng
    integer                          :: indQ, iAt, iL, iM, indLM
    real(dp), allocatable            :: qTrans(:)
    real(dp), allocatable            :: vTmp(:)
    real(dp), allocatable            :: xcGammaU(:,:)
    real(dp), allocatable            :: eeGammaU(:,:)
    real(dp), allocatable            :: tmpU(:,:)
    real(dp), allocatable            :: delQ(:)
    real(dp)                         :: rTmp

    ASSERT(size(specie) == nAtom)
    ASSERT(all(shape(coord) == (/ 3, nAtom /)))
    ASSERT(all(shape(XCSqr) == (/ nOrb, nOrb /)))
    ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    nType = maxval(specie)
    maxM = (maxval(mAngAtom)+1)**2
    ASSERT(all(shape(eTab) == (/ nType, maxM, maxM /)))
    ASSERT(all(shape(hubbU) == (/ maxval(mAngAtom+1), nType /)))
    ASSERT(all(shape(qOutput) == (/ maxval(mAngAtom+1)**2, nAtom /)))
    ASSERT(all(shape(referenceN0) == (/ maxval(mAngAtom+1), nType /)))
    ASSERT(size(eXC) == nOrb)
    
    ALLOCATE_(qTrans, (nAng))
    ALLOCATE_(vTmp, (nOrb))
    ALLOCATE_(xcGammaU, (nAng,nAng))
    ALLOCATE_(eeGammaU, (nAng,nAng))
    ALLOCATE_(tmpU, (maxval(mAngAtom+1),nType))
    ALLOCATE_(delQ, (nAng))

    !! Non SCC part:  <psi_i|v_xc(rho_0)|psi_i>
    eXC(:) = 0.0_dp
    do iOrb = 1, nOrb
      call symv(vTmp, XCSqr, eigenVec(:,iOrb), 'L', 1.0_dp, 0.0_dp)
      eXC(iOrb) = eXC(iOrb) + dot_product(eigenVec(:,iOrb), vTmp)
    enddo
    
    !! Calculate mean electron repulsion integral, equivalent to
    !! atom-resolved Gamma.
    !! Needs to be updated for l-resolved Gamma in SCC
    do iType = 1, nType
      rTmp = 0.0_dp
      mmAng = (mAngSpecie(iType)+1)**2
      do iM1 = 1, mmAng
        do iM2 = iM1, mmAng
          rTmp = rTmp + eTab(iType, iM1, iM2)
        enddo
      enddo
      rTmp = 2.0_dp*rTmp/(mmAng*(mmAng + 1))
      tmpU(:, iType) = rTmp
    enddo
    call gammaM(nAtom, mAngAtom, specie, coord, tmpU, eeGammaU)
    do iType = 1, nType
      tmpU(:,iType) = hubbU(1,iType)
    enddo
    call gammaM(nAtom, mAngAtom, specie, coord, tmpU, xcGammaU)
    xcGammaU(:,:) = xcGammaU(:,:) - eeGammaU(:,:)
    
    !! Charge per L
    delQ(:) = 0.0_dp
    indQ = 0
    do iAt = 1, nAtom
      indLM = 0
      do iL = 1, mAngAtom(iAt)+1
        indQ = indQ + 1
        do iM = 1, 2*iL -1 
          indLM = indLM + 1
          delQ(indQ) = delQ(indQ) + qOutput(indLM, iAt)
        enddo
        delQ(indQ) = delQ(indQ) - referenceN0(iL, specie(iAt))
      enddo
    enddo
    
    !! SCC part: q^i_m * [v_xc]_mn * dq_n; v_xc = Gamma[U_H] - Gamma[U_ee]
    do iOrb = 1, nOrb
      call transQ(nAtom, mAngAtom, iOrb, iOrb, eigenVec, qTrans)
      call symv(vTmp(1:nAng), xcGammaU, delQ, 'L', 1.0_dp, 0.0_dp)
      eXC(iOrb) = eXC(iOrb) + dot_product(qTrans, vTmp(1:nAng))      
    enddo
    
    DEALLOCATE_(delQ)
    DEALLOCATE_(tmpU)
    DEALLOCATE_(vTmp)
    DEALLOCATE_(xcGammaU)
    DEALLOCATE_(eeGammaU)
    DEALLOCATE_(qTrans)
  end subroutine dftbXC

end module gwXandXC
