module GWDriver
#include "allocate.h"  
#include "assert.h"
  use GWGammo
  use GWEpsInv
  use GWPPMReal
  use GWEhfx
  use GWSelfCor
  implicit none
  private

  public :: calculateGW

  
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


  subroutine calculateGW(nAtom, nOrb, mAngAtom, mAngSpecie, mmAng, HSqrReal, &
      &eigen, SSqrReal, &
      &filling, coord0, latVecs, tPeriodic, HSqrReal2, specie0, nSpecie, nEl, &
      &q0Specie, hubbU, xtab, etab, ceri, mass, iAtomStart, qq, ehf, exc, ezse)
    integer, intent(in) :: nAtom
    integer, intent(in) :: nOrb
    integer, intent(in) :: mAngAtom(:)
    integer, intent(in) :: mAngSpecie(:)
    integer, intent(in) :: mmAng
    real(8), intent(in) :: HSqrReal(:,:)
    real(8), intent(in) :: eigen(:)
    real(8), intent(in) :: SSqrReal(:,:)
    real(8), intent(in) :: filling(:)
    real(8), intent(in) :: coord0(:,:)
    real(8), intent(in) :: latVecs(:,:)
    logical, intent(in) :: tPeriodic
    real(8), intent(in) :: HSqrReal2(:,:)
    integer, intent(in) :: specie0(:)
    integer, intent(in) :: nSpecie
    real(8), intent(in) :: nEl
    real(8), intent(in) :: q0Specie(:)
    real(8), intent(in) :: hubbU(:)
    real(8), intent(in) :: xtab(:,:,:)
    real(8), intent(in) :: etab(:,:,:)    
    real(8), intent(in) :: ceri(:,:)
    real(8), intent(in) :: mass(:)
    integer, intent(in) :: iAtomStart(:)
    real*8, intent(in) :: qq(:)
    real*8, intent(out) :: ehf(:)
    real*8, intent(out) :: exc(:)
    real*8, intent(out) :: ezse(:)

    integer, allocatable :: ind(:), lind(:)
    integer :: iAt, nAng
    real(8) :: dacc, basis(3, 3)

    ALLOCATE_(ind, (nAtom+1))
    ind(1:nAtom) = iAtomStart(1:nAtom) - 1
    ind(nAtom+1) = ind(nAtom) + (mAngAtom(nAtom)+1)**2
    ALLOCATE_(lind, (nAtom+1))
    lind(1) = 0
    do iAt = 1, nAtom
      lind(iAt+1) = lind(iAt) + mAngAtom(iAt) + 1
    end do
    nAng = lind(nAtom+1)

    dacc = 4.0_8 * epsilon(1.0_8)
    basis = reshape(latVecs, (/3, 3/), order=(/2, 1/))
    call gw_drv(nAtom, nOrb, (mmAng+1)**2, mmAng+1, ind, HSqrReal, eigen, SSqrReal, &
        & filling, &
        &dacc, coord0, basis, tPeriodic, HSqrReal2, specie0, mAngSpecie+1, &
        &nSpecie, &
        &nEl, q0Specie, hubbU, xtab, etab, ceri, mass, nAng, lind, qq, &
        &ehf, exc, ezse)

    DEALLOCATE_(ind)
    DEALLOCATE_(lind)
    
  end subroutine calculateGW
  



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


  subroutine gw_drv(nn,ndim,morb,ldim, ind,c,ev,s,occ,dacc,rat,basis,period,hxc,izp,&
      &lmax,ntype,rnel,qzero,uhubb,xtab,etab,ceri,xm,nang,lind, qq, &
      &ehf, exc, ezse)
    integer, intent(in) :: nn, ndim, morb, ldim ! morb = mmAng+1
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
    real*8, intent(in) :: qq(:)
    real*8, intent(out) :: ehf(:)
    real*8, intent(out) :: exc(:)
    real*8, intent(out) :: ezse(:)

    !     Local variables
    integer nel,nhel,lind(nn+1),i,j,k,l,m,n,izpj
    real*8 wij(ndim,ndim),epsi(nAng,nAng),zq(nAng),qij(nAng)
    real*8 esc(ndim)
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
    ASSERT(all(shape(xtab) == (/ ntype, morb, morb /)))
    ASSERT(all(shape(etab) == (/ ntype, morb, morb /)))
    ASSERT(all(shape(ceri) == (/ ntype, ldim /)))
    ASSERT(size(xm) == ntype)
    ASSERT(size(ehf) == ndim)
    ASSERT(size(exc) == ndim)
    ASSERT(size(ezse) == ndim)

    write(*,*) '===> GW part starts'
    nel=int(rnel)
    if(mod(nel,2).ne.0) then
      write(*,*) 'Warning: Not a closed shell system!'
    endif
    nhel = nel/2
    if(ev(nhel+1)-ev(nhel).lt.degtol) then
      write(*,*) 'Detected fractional occupation!'
    endif



    !     s-p energy differences
    do i = 1,ndim
      do j = i+1,ndim
        wij(i,j) = ev(j)-ev(i)
      enddo
    enddo

    !     overlap*coefficients needed in routine transq
    call dsymm('L','L',ndim,ndim,1.0d0,s,ndim,c,ndim,0.0d0,stc,ndim)

    !     Get matrix of the coulomb interaction [v] (equal to gamma)
    !     Warning: this gamma is angular momentum resolved 
    call gammo(nn,rat,izp,basis,period,ceri,lind,lmax,gmo,nType,nAng)
    do n=1,nang
      do m=n+1,nang
        gmo(n,m) = gmo(m,n)
      enddo
    enddo


    write(*,*) '===> Epsinv 1. frequency' 
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

    write(*,*) '===> Epsinv 2. frequency' 
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

    write(*,*) '===> X + XC'
    !     Get hartree exchange (ehf) and DFT exchange-correlation (exc)
    call ehfx(nn,ndim,nang,ind,izp,rat,basis,period, &
        &  c,hxc,stc,occ,dacc,qij,gmo,ehf,exc, &
    &  qzero,uhubb,xtab,etab,xm,ntype,lmax,morb, qq)


    write(*,*) '===> Selfcor'
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
    write(*,*) '===> Selfcor epsilon'
    call selfcor(nn,ndim,nang,stc,c,occ,dacc, &
        &  qij,phi,zq,wq,ev,esc,slfc,2,lmax,izp,nType)
    do i = 1,ndim
      sreno(i) = (sreno(i) -  slfc(i))/rdse
      sreno(i) = 1.d0/(1.d0-sreno(i))
      qpcor(i) = sreno(i)*(slfc(i) +ehf(i)-exc(i))
      ezse(i)  = ev(i)+dreal(qpcor(i))
    enddo
    
  end subroutine gw_drv

end module GWDriver
