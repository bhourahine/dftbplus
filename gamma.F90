!=============================================================================
!
!      gamma resulting from exact Coulomb interaction of normalized
!      exp(-a*r) charge distribution
!      Attention: subroutine gamsub needed
!
!      input:  r:     distance
!         uhub1: Hubbard parameter orbital 1
!         uhub2: Hubbard parameter orbital 2
!      output: gval:  gamma12 function value
!
!=============================================================================

       subroutine gam12(r,uhub1,uhub2,gval)
       IMPLICIT NONE
       REAL*8 zero
       parameter(zero=1.0d-4)
       REAL*8 gval,a1,a2,src,avg,uhub1,uhub2,rrc,rrc3
       REAL*8 val12,val21,drv12,drv21,r,fac,fac2,efac
       gval= 0.0
       a1= 3.2*uhub1
       a2= 3.2*uhub2
       IF (a1+a2 .lt. zero) THEN
        RETURN
       ENDIF
       src= 1.0/(a1+a2)
       fac= a1*a2*src
       avg= 1.6*(fac+fac*fac*src)
       IF (r .lt. zero) THEN
           gval= 0.3125*avg
       ELSE
         rrc= 1.0/r
         rrc3= rrc*rrc*rrc
         IF (abs(a1-a2) .lt. 1.0d-5) THEN
           fac= avg*r
           fac2= fac*fac
           efac= exp(-fac)/48.0
           gval= (1.0-(48.0+33*fac+fac2*(9.0+fac))*efac)*rrc
         ELSE
           call gamsub(a1,a2,r,rrc,val12,drv12)
           call gamsub(a2,a1,r,rrc,val21,drv21)
           gval= rrc-val12-val21
         ENDIF
       ENDIF
       RETURN
       END



!=============================================================================
!
!      derivative of gamma resulting from exact Coulomb interaction 
!      of normalized  exp(-a*r) charge distribution
!      Attention: subroutine gamsub needed
!
!      input:  r:     distance
!         uhub1: Hubbard parameter orbital 1
!         uhub2: Hubbard parameter orbital 2
!      output: gdrv:  ((d gamma12)/(d r)) * 1/r
!
!=============================================================================

       subroutine gam121(r,uhub1,uhub2,gdrv)
       IMPLICIT NONE
       REAL*8 zero
       parameter(zero=1.0d-4)
       REAL*8 gdrv,a1,a2,src,avg,uhub1,uhub2,rrc,rrc3
       REAL*8 val12,val21,drv12,drv21,r,fac,fac2,efac
       gdrv= 0.0
       a1= 3.2*uhub1
       a2= 3.2*uhub2
       IF (a1+a2 .lt. zero) THEN
        RETURN
       ENDIF
       src= 1.0/(a1+a2)
       fac= a1*a2*src
       avg= 1.6*(fac+fac*fac*src)
       IF (r .lt. zero) THEN
       ELSE
         rrc= 1.0/r
         rrc3= rrc*rrc*rrc
         IF (abs(a1-a2) .lt. 1.0d-5) THEN
           fac= avg*r
           fac2= fac*fac
           efac= exp(-fac)/48.0
           gdrv= -(1.0-(48.0+48*fac+fac2*(24.0+7*fac+fac2))*efac)*rrc3
         ELSE
           call gamsub(a1,a2,r,rrc,val12,drv12)
           call gamsub(a2,a1,r,rrc,val21,drv21)
           gdrv= -rrc3-(drv12+drv21)*rrc
         ENDIF
       ENDIF
       RETURN
       END



!===========================================================================
!
!      auxiliary routine needed by gam12 and gam121
!
!      input   a:    alpha1
!              b:    alpha2
!              r:    distance
!              rrc:  1/distance
!      output: gval: function value
!              gdrv: function derivative
!
!===========================================================================

       subroutine gamsub(a,b,r,rrc,gval,gdrv)
       IMPLICIT NONE
       REAL*8 a,a2,b,b2,b4,b6,drc,drc2,r,efac,rrc,fac,gval,gdrv
       a2= a*a
       b2= b*b
       b4= b2*b2
       b6= b4*b2
       drc= 1.0/(a2-b2)
       drc2=drc*drc
       efac= exp(-a*r)
       fac= (b6-3*a2*b4)*drc2*drc*rrc
       gval= efac*(0.5*a*b4*drc2-fac)
       gdrv= -a*gval+efac*fac*rrc
       RETURN
       END
