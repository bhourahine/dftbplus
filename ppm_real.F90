!===============================================================================
!
!    ppm_real implements a global diagonalizable plasmon-pole model
!    for [W-v]. This quantity is symmetric in contrast to [eps]. The two
!    test frequencies of the model are pure imaginary.
! 
!    Ansatz: [W-v](w) = phi*L(w)*phi^T; with the eigenvectors phi of [W-v](w=0)
!            and L(w) = zq*wq^2/(w^2-wq^2)
! 
!    Input:
!           wmv    [W-v](w)
!           nang   dimension of GW matrices
!           omega  w
!           mode   1: [W-v] is diagonalized, zq determined
!                  2: wq is determined
!    Output:
!           phi    see above
!           zq     see above
!           wq     see above 
!
!===============================================================================
      subroutine ppm_real(wmv,phi,zq,wq,nang,omega,mode)

      implicit none
      include 'maxima.inc'

!     Parameter
      integer    nang,mode
      real*8     wmv(ADIM,*),phi(ADIM,*),zq(*)
      complex*16 wq(*),omega

!     Local
      integer    i,m,n,info
      real*8     wrk(3*ADIM),ev(ADIM),iomega,tmp

      if(mode.eq.1) then
         if(abs(omega).gt.1.d-10) then
            print *,'ppm_real: omega must be zero for mode = 1'
         endif
         CALL DSYEV('V','U',nang,wmv,ADIM,ev,wrk,3*ADIM,info)
         do m = 1,nang
            do n = 1,nang
               phi(n,m) = wmv(n,m)
            enddo
         enddo 
         do i = 1,nang
            zq(i) = -ev(i)
         enddo
      else
         if(dreal(omega).gt.1.d-10) then
           print *,'ppm_real: omega must be pure imaginary for mode = 2'
         endif
         iomega = dimag(omega)
         CALL DSYEV('N','U',nang,wmv,ADIM,ev,wrk,3*ADIM,info)
         do i = 1,nang
            tmp = -ev(i)*iomega**2/(zq(i)+ev(i))
            if(tmp.lt.0.d0) then
               wq(i) = (0.d0,1.d0)*sqrt(-tmp)
            else
               wq(i) = (1.d0,0.d0)*sqrt(tmp)
            endif
         enddo
      endif
      return
      end




