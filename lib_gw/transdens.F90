!!* Contains subroutines for the calculation of Mulliken populations
!!* representing KS transition densities 
module TransDens
# include "assert.h"
# include "allocate.h"
  use accuracy
  use blasRoutines
  implicit none

  private
  public transQ_init, transQ, transQ_destr
  
  interface transQ_init
    module procedure transQ_init_real
!!    module procedure TransQ_kpts
  end interface
  
  interface transQ
    module procedure transQ_real
!!    module procedure TransQ_kpts
  end interface

  !! Private module variables (suffixed with "_" for clarity)
  logical               :: tTransQinit_ = .false.
  real(dp), allocatable :: sTimEvec_(:,:) ! Overlap times eigenvectors
  
contains

  !!* Initializes transq by evaluating S*C (real case, square multiply)
  !!* @param SSqrReal(nOrb,nOrb)  Overlap matrix (square format)
  !!* @param eigenVec(nOrb,nOrb)  Eigenvectors of Hamiltonian
  subroutine transQ_init_real(SSqrReal,eigenVec)
    real(dp),            intent(in)  :: SSqrReal(:,:)
    real(dp),            intent(in)  :: eigenVec(:,:)
    ASSERT(.not. tTransQinit_)
    ASSERT(all(shape(SSqrReal) == shape(eigenVec)))

    ALLOCATE_(sTimEvec_, (size(SSqrReal, dim=1), size(SSqrReal, dim=2)))
    call symm(sTimEvec_,'L',SSqrReal,eigenVec,'L')
    tTransQinit_ = .true.
  end subroutine transQ_init_real

  !!* Computes l-shell resolved transition charges between Orbital 
  !!* iOrb and jOrb (real case)
  !!* @param nAtom                Number of Atoms
  !!* @param mAngAtom(nAtom)      Max. angular momentum per atom
  !!* @param iOrb                 Kohn-Sham orbital index
  !!* @param jOrb                 Kohn-Sham orbital index
  !!* @param eigenVec(nOrb,nOrb)  Eigenvectors of Hamiltonian
  !!* @param qTrans(nsumL)        Transition charges, l-shell resolved  
  subroutine transQ_real(nAtom,mAngAtom,iOrb,jOrb,eigenVec,qTrans)
    integer,             intent(in)  :: nAtom
    integer,             intent(in)  :: mAngAtom(:)
    integer,             intent(in)  :: iOrb
    integer,             intent(in)  :: jOrb
    real(dp),            intent(in)  :: eigenVec(:,:)
    real(dp),            intent(out) :: qTrans(:)

    integer                          :: nOrb
    integer                          :: iAt,iL,iM,indE,indQ
    
    ASSERT(tTransQinit_)
    ASSERT(size(mAngAtom) == nAtom)
    ASSERT(size(eigenVec,dim=1) == size(eigenVec,dim=2)) 
    nOrb=size(eigenVec,dim=1)
    ASSERT(iOrb <= nOrb)
    ASSERT(jOrb <= nOrb)

    indQ = 0
    indE = 0
    do iAt = 1,nAtom
      do iL = 0,mAngAtom(iAt)
        indQ = indQ + 1
        qTrans(indQ) = 0.0_dp
        do iM = 0,2*iL
          indE = indE + 1
          qTrans(indQ) = qTrans(indQ) + &
              &  0.5_dp*(eigenVec(indE,iOrb)*sTimEvec_(indE,jOrb) &
              &       +  eigenVec(indE,jOrb)*sTimEvec_(indE,iOrb))
        enddo
      enddo
    enddo
  end subroutine transQ_real
    
  !!* Deallocates temp array
  subroutine transQ_destr
    DEALLOCATE_(sTimEvec_)
    tTransQinit_ = .false.
  end subroutine transQ_destr
  
end module TransDens
