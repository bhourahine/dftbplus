!===============================================================================
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
!===============================================================================

      subroutine gw_drv(nn,ndim,ind,c,ev,s,occ,dacc,rat,basis,period, &
                     &  hxc,izp,lmax,rnel,ceri,ntype,qzero,uhubb,xtab, &
                     &  etab,xm,lpol) 

      implicit none
      include 'maxima.inc'

!     Conversion hartree -> eV 
      real*8 conv
      parameter (conv = 27.21139908d0)
!     Step for differentiation of selfenergy
      real*8 rdse
      parameter (rdse = 5.d-4)
!     degeneracy tolerance for T = 0
      real*8 degtol 
      parameter (degtol = 1.0d-4) 

!     Parameter
      integer nn,ndim,ind(*)
      real*8  c(MDIM,*),ev(MDIM),s(MDIM,*),occ(*),dacc
      real*8  rat(3,*),basis(3,*),hxc(MDIM,*)
      logical period

!     Local variables
      integer nel,nhel,lind(NNDIM+1),nang,i,j,k,l,m,n,izpj
      real*8 wij(MDIM,MDIM),epsi(ADIM,ADIM),zq(ADIM),qij(ADIM)
      real*8 esc(MDIM),ehf(MDIM),exc(MDIM),ezse(MDIM)
      real*8 wmv(ADIM,ADIM),phi(ADIM,ADIM),gmo(ADIM,ADIM)
      real*8 stc(MDIM,MDIM),wzero(ADIM,ADIM)
      complex*16 omega,wq(ADIM),slfc(MDIM),sreno(MDIM),qpcor(MDIM)

!     Common block variables 
      integer izp(NNDIM),nbeweg
      integer lmax(MAXTYP),ntype
      real*8 rnel,qzero(MAXTYP)
      real*8 uhubb(MAXTYP),uhbuu(MAXTYP),uhbud(MAXTYP)
      real*8 xtab(MAXTYP,LDIM,LDIM),etab(MAXTYP,LDIM,LDIM)
      real*8 ceri(MAXTYP,LDIM),xm(MAXTYP)
      logical lpol


      print *,'===> GW part starts'
      nel=int(rnel)
      if(mod(nel,2).ne.0) then
        print *,'Warning: Not a closed shell system!'
      endif
      nhel = nel/2
      if(ev(nhel+1)-ev(nhel).lt.degtol) then
        print *,'Detected fractional occupation!'
      endif


!     Dimension of GW matrices
      lind(1) = 0
      do j = 1,nn 
         izpj = izp(j)
         lind(j+1) = lind(j)+lmax(izpj)
      end do
      nang = lind(nn+1)
      
!     s-p energy differences
      do i = 1,ndim
        do j = i+1,ndim
          wij(i,j) = ev(j)-ev(i)
        enddo
      enddo

!     overlap*coefficients needed in routine transq
      call dsymm('L','U',ndim,ndim,1.0d0,s,MDIM,c,MDIM,0.0d0,stc,MDIM)

!     Get matrix of the coulomb interaction [v] (equal to gamma)
!     Warning: this gamma is angular momentum resolved 
      call gammo(nn,rat,izp,basis,period,ceri,lind,lmax,gmo)
      do n=1,nang
        do m=n+1,nang
           gmo(n,m) = gmo(m,n)
        enddo
      enddo
      

      print *,'===> Epsinv 1. frequency' 
!     Set up plasmon-pole model, epsi = <eps>^-1, wmv = epsi*v - v  
      omega = (0.d0,0.d0)     
      call epsinv(nn,ndim,nang,gmo,stc,c,occ,dacc, &
                  &  qij,wij,omega,epsi,lmax,izp)
      do m = 1,nang
         do n = m,nang
            wmv(n,m) = gmo(n,m)
            wmv(m,n) = gmo(n,m)
         enddo
      enddo
      CALL DSYMM('R','L',nang,nang,1.d0,gmo,ADIM, &
                      &  epsi,ADIM,-1.d0,wmv,ADIM)
!     Calculate W for zero frequency
      CALL DSYMM('R','L',nang,nang,1.d0,gmo,ADIM, &
                      &  epsi,ADIM,0.d0,wzero,ADIM)

      call ppm_real(wmv,phi,zq,wq,nang,omega,1)

      print *,'===> Epsinv 2. frequency' 
!     Do same thing for another testfrequecy
      omega = (0.d0,5.d-1)
      call epsinv(nn,ndim,nang,gmo,stc,c,occ,dacc, &
                 &  qij,wij,omega,epsi,lmax,izp)
      do m = 1,nang
         do n = m,nang
            wmv(n,m) = gmo(n,m)
            wmv(m,n) = gmo(n,m)
         enddo
      enddo
      CALL DSYMM('R','L',nang,nang,1.d0,gmo,ADIM, &
                      &  epsi,ADIM,-1.d0,wmv,ADIM)
      call ppm_real(wmv,phi,zq,wq,nang,omega,2)
!     PPM parameters w_q and z_q are now known, phi are eigenvectors of W-v

      print *,'===> X + XC'
!     Get hartree exchange (ehf) and DFT exchange-correlation (exc)
      call ehfx(nn,ndim,nang,ind,izp,rat,basis,period, &
             &  c,hxc,stc,occ,dacc,qij,gmo,ehf,exc,
             &  qzero,uhubb,xtab,etab,xm,ntype,lmax)


      print *,'===> Selfcor'
!     Consider energy dependence of selfenergy through renormalization
!     Selfenergy is numerically differentiated 
      do i = 1,ndim
         esc(i) = ev(i) + rdse
      enddo
      call selfcor(nn,ndim,nang,stc,c,occ,dacc, &
                  &  qij,phi,zq,wq,ev,esc,slfc,1,lmax,izp)
      do i = 1,ndim
         sreno(i) = slfc(i)
      enddo
      call DCOPY(ndim,ev,1,esc,1)
      print *,'===> Selfcor epsilon'
      call selfcor(nn,ndim,nang,stc,c,occ,dacc, &
                  &  qij,phi,zq,wq,ev,esc,slfc,2,lmax,izp)
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
      end
