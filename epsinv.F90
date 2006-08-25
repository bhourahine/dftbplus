!===============================================================================
!
!    epsinv calculates ([eps(w)]S)^-1 = (1 - [v]<P(w)>)^-1
!   
!    Input: 
!           nn     number of atoms
!           ndim   number of basis functions
!           nang   dimension of GW matrices
!           gam    [v]
!           stc    overlap matrix times coefficients
!           c      MO coeffiecients
!           occ    occupation numbers
!           dacc   occupation number tolerance
!           qij    transition charges between orbitals i and j
!           wij    MO energy differences between orbitals i and j
!           omega  frequency w
!    Output:
!           epsi    ([eps(w)]S)^-1 
!
!===============================================================================
      subroutine epsinv(nn,ndim,nang,gam,stc,c,occ,dacc, &
                      & qij,wij,omega,epsi,lmax,izp)

      implicit none
      include 'maxima.inc'

      real*8 epscor
      parameter (epscor = 1.d0)
!     broadening for polarization
      real*8 brd
      parameter (brd = 1.d-8)

!     Parameter
      integer    nn,ndim,nang
      real*8     gam(ADIM,*),stc(MDIM,*),c(MDIM,*),occ(*),dacc 
      real*8     qij(*),wij(MDIM,*),epsi(ADIM,*)
      complex*16 omega
      integer  izp(NNDIM),lmax(MAXTYP)

!     Local
      integer    n,m,info,ipiv(ADIM)
      real*8     wrk(64*ADIM),plr(ADIM,ADIM),epsmac,tmp



      if(abs(dreal(omega)).gt.1.d-10) then
         print *,'EPSINV works only for imaginary omega'
         stop
      endif

!     Build polarization matrix <P(w)> 
!     We only need real part of <P>, because we use
!     pure imaginary frequencies to build up the plasmon pole model
      call polar(nn,ndim,nang,stc,c,occ,dacc,qij,wij,omega,plr,brd,lmax,izp)
      do m = 1,nang
         do n = 1,m-1
            epsi(m,n) = 0.d0
            epsi(n,m) = 0.d0
         enddo
         epsi(m,m) = 1.d0
      enddo

      CALL DSYMM('L','L',nang,nang,-epscor,gam, &
              &   ADIM,plr,ADIM,1.d0,epsi,ADIM)

      tmp = 0.d0
      do m = 1,nang
         tmp = tmp + epsi(m,m)
      enddo
      write(*,'(x,A,f10.6)') 'Eps: ',tmp/nang

      CALL DGETRF(nang,nang,epsi,ADIM,ipiv,info)
      CALL DGETRI(nang,epsi,ADIM,ipiv,wrk,64*ADIM,info)

      return
      end









