!!* Contains subroutines to calculate the RPA polarizabilty,
!!* dielectric constant and inverse of it, needed in GW calculation
module RandomPhase
# include "assert.h"
# include "allocate.h"
  use accuracy, only : dp, brdRPA
  use blasRoutines
  use lapackRoutines
  use TransDens
  implicit none

  private
  public buildPolar, buildEpsInv
  
contains

  !!* Computes the real part of the RPA polarizability <P(w)> for purely
  !!* imaginary w
  !!* @param nAtom                  Number of Atoms
  !!* @param nOrb                   Number of orbitals
  !!* @param nAng                   Size of GW basis
  !!* @param mAngAtom(nAtom)        Max. angular momentum per atom
  !!* @param filling(nOrb)          Occupation of MO's (spin restricted)
  !!* @param eigenVec(nOrb,nOrb)    Eigenvectors of Hamiltonian
  !!* @param delEigenVal(nOrb,nOrb) Eigenvalue differences
  !!* @param omegaPPi               Imaginary part of w
  !!* @param polar(nAng,nAng)       <P(w)> polarizability, symmetric 
  subroutine buildPolar(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec,&
      & delEigenVal, omegaPPi, polar)
    integer,             intent(in)  :: nAtom
    integer,             intent(in)  :: nOrb
    integer,             intent(in)  :: nAng
    integer,             intent(in)  :: mAngAtom(:)
    real(dp),            intent(in)  :: filling(:)
    real(dp),            intent(in)  :: eigenVec(:,:)
    real(dp),            intent(in)  :: delEigenVal(:,:)
    real(dp),            intent(in)  :: omegaPPi
    real(dp),            intent(out) :: polar(:,:)

    integer                          :: iOrb, jOrb, iAng, jAng
    real(dp)                         :: delOcc
    real(dp), allocatable            :: qTrans(:)

    ASSERT(size(mAngAtom) == nAtom)
    ASSERT(size(filling) == nOrb)
    ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    ASSERT(all(shape(delEigenVal) == (/ nOrb, nOrb /)))
    ASSERT(all(shape(polar) == (/ nAng, nAng /)))
    ALLOCATE_(qTrans, (nAng))
    polar(:,:) = 0.0_dp
     
    do iOrb = 1, nOrb
      do jOrb = iOrb+1, nOrb
        call transQ(nAtom, mAngAtom, iOrb, jOrb, eigenVec, qTrans)
        delOcc = filling(iOrb) - filling(jOrb)
        if(abs(delOcc) >  epsilon(1.0_dp)) then
          do iAng = 1, nAng
            do jAng = iAng, nAng
              polar(jAng,iAng) = polar(jAng,iAng) &
                  &-2.0_dp*delOcc*qTrans(jAng)*qTrans(iAng)*&
                  &(delEigenVal(iOrb,jOrb) + brdRPA)/&
                  &(omegaPPi**2 + (delEigenVal(iOrb,jOrb) + brdRPA)**2)
            enddo
          enddo
        endif
      enddo
    enddo
    do iAng = 1, nAng
      do jAng = iAng+1, nAng
        polar(iAng,jAng) =  polar(jAng,iAng)
      enddo
    enddo
    
    DEALLOCATE_(qTrans)
  end subroutine buildPolar
  
  !!* Computes RPA dielectric matrix and inverse of it:
  !!* ([eps(w)]S)^-1 = (1 - [v]<P(w)>)^-1
  !!* @param nAtom                  Number of Atoms
  !!* @param nOrb                   Number of orbitals
  !!* @param nAng                   Size of GW basis
  !!* @param mAngAtom(nAtom)        Max. angular momentum per atom
  !!* @param filling(nOrb)          Occupation of MO's (spin restricted)
  !!* @param eigenVec(nOrb,nOrb)    Eigenvectors of Hamiltonian
  !!* @param delEigenVal(nOrb,nOrb) Eigenvalue differences
  !!* @param eeGamma(nAng,nAng)     [v]
  !!* @param omegaPPi               Imaginary part of w
  !!* @param epsInv(nAng,nAng)      ([eps(w)]S)^-1
  subroutine buildEpsInv(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec,&
      & delEigenVal, eeGamma, omegaPPi, epsInv)
    integer,             intent(in)  :: nAtom
    integer,             intent(in)  :: nOrb
    integer,             intent(in)  :: nAng
    integer,             intent(in)  :: mAngAtom(:)
    real(dp),            intent(in)  :: filling(:)
    real(dp),            intent(in)  :: eigenVec(:,:)
    real(dp),            intent(in)  :: delEigenVal(:,:)
    real(dp),            intent(in)  :: eeGamma(:,:)
    real(dp),            intent(in)  :: omegaPPi
    real(dp),            intent(out) :: epsInv(:,:)

    integer                          :: iAng
    integer, allocatable             :: ipiv(:)
    real(dp), allocatable            :: polar(:,:)
    
    
    ASSERT(size(mAngAtom) == nAtom)
    ASSERT(size(filling) == nOrb)
    ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    ASSERT(all(shape(delEigenVal) == (/ nOrb, nOrb /)))
    ASSERT(all(shape(eeGamma) == (/ nAng, nAng /)))
    ASSERT(all(shape(epsInv) == (/ nAng, nAng /)))
    ALLOCATE_(polar, (nAng,nAng))
    ALLOCATE_(ipiv, (nAng))

    call buildPolar(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec,&
        & delEigenVal, omegaPPi, polar)
    epsInv(:,:) = 0.0_dp
    do iAng = 1, nAng
      epsInv(iAng,iAng) = 1.0_dp
    enddo
    ! Compute dielectric matrix: 1 - [v]<P(w)>
    call symm(epsInv, 'L', eeGamma, polar, 'L', -1.0_dp, 1.0_dp)
    ! Taking the inverse
    call getrf(epsInv, ipiv)
    call getri(epsInv, ipiv)
    
    DEALLOCATE_(ipiv)   
    DEALLOCATE_(polar)
  end subroutine buildEpsInv
end module RandomPhase
