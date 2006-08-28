module GWSelfCor
  use GWTransQ
  implicit none
  private

  public :: selfcor

contains
  
  
  !===============================================================================
  !
  !    selfcor calculates the correlation part of the self energy 
  !    <psi_i|Sigma_c(E)|psi_i>
  !
  !    Input: 
  !           nn     number of atoms
  !           ndim   number of basis functions
  !           nang   dimension of GW matrices 
  !           occ    occupation numbers
  !           dacc   occupation number tolerance
  !           stc    overlap matrix times coeffcients
  !           c      MO coeffiecients
  !           qij    transition charges between orbitals i and j
  !           phi    eigenvectors phi of [W-v](w=0)
  !           zq     parameters of the plasmon pole model (see ppm_real)
  !           wq     parameters of the plasmon pole model (see ppm_real)
  !           ev     KS eigenvalues
  !           wij    MO energy differences between orbitals i and j
  !           omega  frequency w
  !           esc    energies at which to evaluate Sigma_c
  !           mode   see below
  !    Output:
  !           qij    qij*phi
  !           slfc   <psi_i|Sigma_c(E)|psi_i> 
  !
  !===============================================================================
  subroutine selfcor(nn,ndim,nang,stc,c,occ,dacc, &
      &  qij,phi,zq,wq,ev,esc,slfc,mode,lmax,izp,nType)

    !     Parameter
    integer    nn,ndim,nang,mode
    real*8     stc(nDim,*),c(nDim,*),occ(*),dacc
    real*8     qij(*),phi(nAng,*),zq(*),ev(*),esc(*)
    complex*16 wq(*),slfc(*)
    integer  izp(nn),lmax(nType)

    !     Local
    integer    i,n,j,k,l,hocc
    real*8     tqij
    integer, intent(in) :: nType

    !     Find highest occupied orbital 
    do k = 1,ndim
      if(abs(occ(k)).lt.dacc) exit
    enddo
    hocc = k-1

    do i = 1,ndim
      slfc(i) = (0.d0,0.d0)
      do n = 1,hocc
        call transq(nn,stc,c,qij,i,n,lmax,izp,ndim,nType)
        do l = 1,nang
          if(zabs(wq(l)).gt.1.d-10) then
            tqij = 0.d0
            do k = 1,nang
              tqij =  tqij + qij(k)*phi(k,l)
            enddo
            slfc(i) =  slfc(i) + 0.5d0*occ(n)*tqij*tqij* &
                & (0.5d0*zq(l)*wq(l)/(esc(i)-ev(n)+wq(l)))
          endif
        enddo
      enddo
      do n = hocc+1,ndim
        call transq(nn,stc,c,qij,i,n,lmax,izp,ndim,nType)
        do l = 1,nang
          if(zabs(wq(l)).gt.1.d-10) then
            tqij = 0.d0
            do k = 1,nang
              tqij =  tqij + qij(k)*phi(k,l)
            enddo
            slfc(i) =  slfc(i) + tqij*tqij* &
                & (0.5d0*zq(l)*wq(l)/(esc(i)-ev(n)-wq(l)))
          endif
        enddo
      enddo
    enddo

  end subroutine selfcor


end module GWSelfCor
