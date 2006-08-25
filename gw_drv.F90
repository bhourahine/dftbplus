module GWDriver
#include "assert.h"
  use GWGammo
  use GWEpsInv
  use GWPPMReal
  use GWEhfx
  use GWSelfCor
  implicit none
  private

  public :: gw_drv

  
  !     Conversion hartree -> eV 
  real*8 conv
  parameter (conv = 27.21139908d0)
  !     Step for differentiation of selfenergy
  real*8 rdse
  parameter (rdse = 5.d-4)
  !     degeneracy tolerance for T = 0
  real*8 degtol 
  parameter (degtol = 1.0d-4) 

  

contains

  !!!     Dimension of GW matrices
  !!  lind(1) = 0
  !!  do j = 1,nn 
  !!    izpj = izp(j)
  !!    lind(j+1) = lind(j)+lmax(izpj)
  !!  end do
  !!  nang = lind(nn+1)


  !=============================================================================
  !
  !    gw_drv calculates quasi particle corrections to the DFTB eigenvalues
  !    using an approximate GW method
  !   
  !    Input: nn     number of atoms
  !           ndim   number of basis functions
  !           ind    indexes AO's belonging to an atom
  !           c      MO coefficients
  !           ev     KS eigenvalues
  !           s      overlap matrix
  !           occ    occupation numbers
  !           dacc   occupation number tolerance
  !           rat    atom positions
  !           basis  basis of cell  
  !           period periodic boundary conditions?
  !           hxc    XC part of the non-SCC hamiltonian
  !
  !============================================================================


  subroutine gw_drv(nn,ndim,ldim,ind,c,ev,s,occ,dacc,rat,basis,period,hxc,izp,&
      &lmax,ntype,rnel,qzero,uhubb,xtab,etab,ceri,xm,nang)
    integer, intent(in) :: nn, ndim, ldim ! ldim = mmAng+1
    integer, intent(in) :: ind(:)
    real*8, intent(in) ::  c(:,:)
    real*8, intent(in) ::  ev(:)
    real*8, intent(in) ::  s(:,:)
    real*8, intent(in) ::  occ(:)
    real*8, intent(in) :: dacc
    real*8, intent(in) :: rat(:,:)
    real*8, intent(in) :: basis(:,:)
    logical, intent(in) :: period
    real*8, intent(in) :: hxc(:,:)
    ! Old common block variables 
    integer, intent(in) ::  izp(:)   ! specie
    integer, intent(in) ::  lmax(:)  ! ang mom per species (s=1)
    integer, intent(in) :: ntype      ! nSpecie
    real*8, intent(in) :: rnel       ! nEl
    real*8, intent(in) :: qzero(:)   ! neutral charges for each species
    real*8, intent(in) :: uhubb(:)   ! Modified Us
    real*8, intent(in) :: xtab(:,:,:)
    real*8, intent(in) :: etab(:,:,:)
    real*8, intent(in) :: ceri(:,:)
    real*8, intent(in) :: xm(:)  ! Mass (for checking for ghosts)
    integer, intent(in) :: nang

    !     Local variables
    integer nel,nhel,lind(nn+1),i,j,k,l,m,n,izpj
    real*8 wij(ndim,ndim),epsi(nAng,nAng),zq(nAng),qij(nAng)
    real*8 esc(ndim),ehf(ndim),exc(ndim),ezse(ndim)
    real*8 wmv(nAng,nAng),phi(nAng,nAng),gmo(nAng,nAng)
    real*8 stc(ndim,ndim),wzero(nAng,nAng)
    complex*16 omega,wq(nAng),slfc(ndim),sreno(ndim),qpcor(ndim)
    
    ASSERT(size(ind) >= nn + 1)
    ASSERT(all(shape(c) == (/ ndim, ndim /)))
    ASSERT(size(ev) == ndim)
    ASSERT(all(shape(s) == (/ ndim, ndim /)))
    ASSERT(size(occ) == ndim)
    ASSERT(all(shape(rat) == (/ 3, nn /)))
    ASSERT(all(shape(basis) == (/3, 3 /)))
    ASSERT(all(shape(hxc) == (/ ndim, ndim /)))
    ASSERT(size(izp) == nn)
    ASSERT(size(lmax) == ntype)
    ASSERT(size(qzero) == ntype)
    ASSERT(size(uhubb) == ntype)
    ASSERT(all(shape(xtab) == (/ nn, ldim, ldim /)))
    ASSERT(all(shape(etab) == (/ nn, ldim, ldim /)))
    ASSERT(all(shape(ceri) == (/ ntype, ldim /)))
    ASSERT(size(xm) == ntype)



    print *,'===> GW part starts'
    nel=int(rnel)
    if(mod(nel,2).ne.0) then
      print *,'Warning: Not a closed shell system!'
    endif
    nhel = nel/2
    if(ev(nhel+1)-ev(nhel).lt.degtol) then
      print *,'Detected fractional occupation!'
    endif



    !     s-p energy differences
    do i = 1,ndim
      do j = i+1,ndim
        wij(i,j) = ev(j)-ev(i)
      enddo
    enddo

    !     overlap*coefficients needed in routine transq
    call dsymm('L','U',ndim,ndim,1.0d0,s,ndim,c,ndim,0.0d0,stc,ndim)

    !     Get matrix of the coulomb interaction [v] (equal to gamma)
    !     Warning: this gamma is angular momentum resolved 
    call gammo(nn,rat,izp,basis,period,ceri,lind,lmax,gmo,nType,nAng)
    do n=1,nang
      do m=n+1,nang
        gmo(n,m) = gmo(m,n)
      enddo
    enddo


    print *,'===> Epsinv 1. frequency' 
    !     Set up plasmon-pole model, epsi = <eps>^-1, wmv = epsi*v - v  
    omega = (0.d0,0.d0)     
    call epsinv(nn,ndim,nang,gmo,stc,c,occ,dacc, &
        &  qij,wij,omega,epsi,lmax,izp,nType)
    do m = 1,nang
      do n = m,nang
        wmv(n,m) = gmo(n,m)
        wmv(m,n) = gmo(n,m)
      enddo
    enddo
    CALL DSYMM('R','L',nang,nang,1.d0,gmo,nAng, &
        &  epsi,nAng,-1.d0,wmv,nAng)
    !     Calculate W for zero frequency
    CALL DSYMM('R','L',nang,nang,1.d0,gmo,nAng, &
        &  epsi,nAng,0.d0,wzero,nAng)

    call ppm_real(wmv,phi,zq,wq,nang,omega,1)

    print *,'===> Epsinv 2. frequency' 
    !     Do same thing for another testfrequecy
    omega = (0.d0,5.d-1)
    call epsinv(nn,ndim,nang,gmo,stc,c,occ,dacc, &
        &  qij,wij,omega,epsi,lmax,izp, nAng)
    do m = 1,nang
      do n = m,nang
        wmv(n,m) = gmo(n,m)
        wmv(m,n) = gmo(n,m)
      enddo
    enddo
    CALL DSYMM('R','L',nang,nang,1.d0,gmo,nAng, &
        &  epsi,nAng,-1.d0,wmv,nAng)
    call ppm_real(wmv,phi,zq,wq,nang,omega,2)
    !     PPM parameters w_q and z_q are now known, phi are eigenvectors of W-v

    print *,'===> X + XC'
    !     Get hartree exchange (ehf) and DFT exchange-correlation (exc)
    call ehfx(nn,ndim,nang,ind,izp,rat,basis,period, &
        &  c,hxc,stc,occ,dacc,qij,gmo,ehf,exc, &
    &  qzero,uhubb,xtab,etab,xm,ntype,lmax,ldim)


    print *,'===> Selfcor'
    !     Consider energy dependence of selfenergy through renormalization
    !     Selfenergy is numerically differentiated 
    do i = 1,ndim
      esc(i) = ev(i) + rdse
    enddo
    call selfcor(nn,ndim,nang,stc,c,occ,dacc, &
        &  qij,phi,zq,wq,ev,esc,slfc,1,lmax,izp,nType)
    do i = 1,ndim
      sreno(i) = slfc(i)
    enddo
    call DCOPY(ndim,ev,1,esc,1)
    print *,'===> Selfcor epsilon'
    call selfcor(nn,ndim,nang,stc,c,occ,dacc, &
        &  qij,phi,zq,wq,ev,esc,slfc,2,lmax,izp,nType)
    do i = 1,ndim
      sreno(i) = (sreno(i) -  slfc(i))/rdse
      sreno(i) = 1.d0/(1.d0-sreno(i))
      qpcor(i) = sreno(i)*(slfc(i) +ehf(i)-exc(i))
      ezse(i)  = ev(i)+dreal(qpcor(i))
    enddo

    open(45,FILE='SEL.DAT')
    write(45,*) & 
        & 'Orb     eps.      E_X         E_XC       E_HF       E_QP'
    write(45,*) &
        & '=============================================================='
    do i = 1,ndim
      write(45,'(1x,i4,3x,5(f8.4,4x))') i,conv*ev(i),conv*ehf(i) &
          &  ,conv*exc(i),conv*(ev(i)-exc(i)+ehf(i)) &
          &  ,conv*(ezse(i))
    enddo
    close(45)
    
    return

  end subroutine gw_drv

end module GWDriver
