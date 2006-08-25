!   ************************************************************************
!   * THIS SUBROUTINE CALCULATES THE RPA POLARIZABILITY                    *
!   *                                                                      *
!   * ON INPUT:                                                            *
!   *           NN:     NUMBER OF ATOMS [I]                                *
!   *           OCC:    OCCUPATION NUMBERS                                 *
!   *           DACC:   OCCUPATION NUMBER TOLERANCE                        *
!   *           QIJ:    TRANSITION CHARGES [R(MDIM,MDIM,*)                 *
!   *           WIJ:    KS-ENERGYDIFFERENCES [R(MDIM,*)]                   *
!   *           OMEGA:  FREQUENCY [C]                                      *
!   *           BRD:    Broadening factor in eV [R]                        *
!   *                                                                      *
!   * ON OUTPUT:                                                           *
!   *           PLR:    POLARIZATION  [R(ADIM,*)]                          *                
!   ************************************************************************

      subroutine polar(nn,ndim,nang,stc,c,occ,dacc, &
                      & qij,wij,omega,plr,brd,lmax,izp)

      implicit none
      include 'maxima.inc'

      integer nn,m,n,ndim,nang,i,j
      real*8 stc(MDIM,*),c(MDIM,*),occ(*),dacc
      real*8 qij(*),wij(MDIM,*),plr(ADIM,*)
      real*8 brd,iomega,focc,tmp
      complex*16 omega
      integer  izp(NNDIM),lmax(MAXTYP)

      if(dreal(omega).gt.1.d-10) then
         print *,'polar: omega must be pure imaginary!'
         stop
      endif
      iomega = dimag(omega)
      do n = 1,nang
         do m = n,nang
            plr(m,n) = 0.d0
         enddo
      enddo
      do i = 1,ndim
         do j = i+1,ndim
            focc = occ(i)-occ(j)
            if(abs(focc).gt.dacc) then
              call transq(nn,stc,c,qij,i,j,lmax,izp)
              do n = 1,nang
                 do m = n,nang
                    plr(m,n) =  plr(m,n) - 2.d0*focc*qij(m)*qij(n) &
                      &  *(wij(i,j)+brd)/(iomega**2+(wij(i,j)+brd)**2)
                 enddo
              enddo
            endif
         enddo
      enddo

      do n = 1,nang
         do m = n+1,nang
            plr(n,m) =  plr(m,n)
         enddo
      enddo
      return
      end 





















