module GWTransQ
  implicit none
  private

  public :: transq

contains
  
  !
  !     computes the transition charges qij between two orbitals i,j
  !
  !     input: 
  !             nn:     number of atoms
  !             ndim:   number of basis functions
  !             ind:    index vector 
  !             stimc:  overlap * MO-coeffcients
  !             c:      MO-coefficients
  !             i,j:    indices of orbitals
  !     output:
  !             qij:    atomic overlap charges
  ! 
  subroutine transq(nn,stimc,c,qij,i,j,lmax,izp,ndim,nType)
    !     Parameter
    integer    nn,i,j
    real*8     stimc(nDim,*),c(nDim,*),qij(*)

    !     Local
    integer    indl,indm,izpn,n,l,m

    !     Common block variables
    integer  izp(nn),lmax(nType)

    integer, intent(in) :: ndim
    integer, intent(in) :: nType


    indl = 0
    indm = 0
    do n = 1,nn
      izpn = izp(n)
      do l = 1,lmax(izpn)
        indl = indl + 1
        qij(indl) = 0.d0  
        do m = 1,2*l-1
          indm = indm + 1
          qij(indl) = qij(indl) + (c(indm,i)*stimc(indm,j) &
              &  +  c(indm,j)*stimc(indm,i))*0.5d0
        enddo
      enddo
    enddo

  end subroutine transq


end module GWTransQ
