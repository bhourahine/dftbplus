!!* Contains subroutines for the calculation of DFTB Gamma Matrices
module GammaMat
# include "assert.h"
# include "allocate.h"
  use accuracy
  use coulomb
  use short_gamma
  implicit none

  private
  public gammaM
  
  interface gammaM
    module procedure gammaM_real
  end interface

contains
  !!* Calculates l-resolved DFTB Gamma matrix (lower triangular) 
  !!* @param nAtom                Number of Atoms
  !!* @param mAngAtom(nAtom)      Max. angular momentum per atom 
  !!* @param specie(nAtom)        species list of atomic species for all atom
  !!* @param coord(3,nAtom)       Atomic coordinates
  !!* @param hubbU(l,nType)       Electron repulsion constants = Hubbard values
  !!* @param gamma(nAng,nAng)     Gamma matrix
  subroutine gammaM_real(nAtom, mAngAtom, specie, coord, hubbU, gamma)
    integer,             intent(in)  :: nAtom
    integer,             intent(in)  :: mAngAtom(:)
    integer,             intent(in)  :: specie(:)
    real(dp),            intent(in)  :: coord(:,:)
    real(dp),            intent(in)  :: hubbU(:,:)
    real(dp),            intent(out) :: gamma(:,:)

    integer                          :: nAng
    integer                          :: iAt1, iAt2, iL1, iL2
    integer                          :: indL(nAtom+1)
    real(dp)                         :: vect(3), r12
    real(dp)                         :: U1, U2
    real(dp)                         :: invRMat(nAtom,nAtom)

    ASSERT(all(shape(coord) == (/ 3, nAtom /)))
    ASSERT(size(hubbU, dim=1) == maxval(mAngAtom+1))
    ASSERT(size(specie) == nAtom)
    nAng = sum(mAngAtom) + nAtom
    ASSERT(all(shape(gamma) == (/ nAng, nAng /)))

    call invR(invRMat, nAtom, coord)
    gamma(:,:) = 0.0_dp
    indL(1) = 0
    do iAt1 = 1, nAtom
      indL(iAt1+1) = indL(iAt1) + mAngAtom(iAt1) + 1
    enddo
    
    do iAt1 = 1, nAtom
      do iL1 = 1, mAngAtom(iAt1)+1
        U1 = hubbU(iL1,specie(iAt1))
        do iAt2 = iAt1, nAtom
          vect(:) = coord(:,iAt1) - coord(:,iAt2)
          r12 = sqrt(sum(vect(:)**2))
          do iL2 = 1, mAngAtom(iAt2)+1
            U2 = hubbU(iL2,specie(iAt2))
            gamma(indL(iAt2)+iL2,indL(iAt1)+iL1) = invRMat(iAt2,iAt1) &
                &                                - expGamma(r12,U2,U1)
          enddo
        enddo
      enddo
    enddo
  end subroutine gammaM_real
  
end module GammaMat
