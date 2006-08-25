!=============================================================================
!   Build lower triangular matrix containing Ewald potentials + short
!   range terms
!
!   INPUT Parameter:
!   INTEGER nat          number of atoms
!   REAL*8 rat(3,*)      position of atoms
!   REAL*8 u(*)          hubbard parameters
!   REAL*8 basis(3,3)    basis of cell
!   LOGICAL period       .true. if periodic boundary conditions
!
!   OUTPUT:
!   REAL*8 gammamat(*,*) matrix containing the values of the ewlad potential
!                      in the upper triangular part
!      !!! NOTE THAT phi(ri - rj) = phi(rj - ri) !!!
!      !!! NOTE THAT shortrange(ri - rj) = shortrange(rj - ri) !!!
!
!=============================================================================

       SUBROUTINE gammamatrix(nat,rat,atomtype,basis,period,u,gammamat)
       IMPLICIT NONE
       INCLUDE 'maxima.inc'
       INTEGER nat
       REAL*8 rat(3,*),basis(3,3),tol,u(MAXTP),gammamat(NNDIM,NNDIM)
       INTEGER i,j,atomtype(*)
       REAL*8 tmpvalue,phivalue,r(3)
       REAL*8 recbasis(3,3), vol
       REAL*8 alpha
       EXTERNAL getalpha
       REAL*8 getalpha
       REAL*8 gval,norm
       EXTERNAL GAM12
       LOGICAL period

       
       do i=1,nat
        do j=1,i
         r(1)=rat(1,i)-rat(1,j)
         r(2)=rat(2,i)-rat(2,j)
         r(3)=rat(3,i)-rat(3,j)

!       !!! NOTE THAT phi(ri - rj) = phi(rj - ri) !!!
!       !!! NOTE THAT shortrange(ri - rj) = shortrange(rj - ri) !!!

       
       IF(period) THEN

!      get reciprocal lattice vectors and cell volume
       CALL REZVOL(basis,recbasis,vol)

!      choose good convergence parameter alpha
       alpha = getalpha(basis)

!      set tolerance for convergence
       tol = 1.0d-8

       CALL phi(r,basis,recbasis,alpha,vol,tol,phivalue) 
       CALL shortrange(r,basis,u(atomtype(i)),u(atomtype(j)), &
                    &  tol,tmpvalue)

       gammamat(i,j)=tmpvalue + phivalue

       ELSE

       norm   = sqrt(r(1)**2+r(2)**2+r(3)**2)

!      get value for Gamma
       CALL GAM12(norm,u(atomtype(i)),u(atomtype(j)),gval)

       gammamat(i,j)=gval

       END IF

       END DO
       END DO
       
       END


!=============================================================================
!   Build lower triang. matrices containg derivatives of long+shortrange expr.
!   INPUT Parameter:
!   INTEGER nat             Number of atoms
!   REAL*8 rat(3,*)         position of atoms
!   REAL*8 u(MAXTP)             hubbard parameters
!   REAL*8 basis(3,3)       basis of cell
!   LOGICAL period          .true. if periodic boundary conditions
!
!   OUTPUT:
!   REAL*8 gammamat1x(*,*)  x derivative of ewald potential in upper
!                           triangular part
!   REAL*8 gammamat1y(*,*)  y derivative of ewald potential in upper
!                           triangular part
!   REAL*8 gammamat1z(*,*)  z derivative of ewald potential in upper
!                           triangular part
!   !!! NOTE THAT shortrange1(ri - rj) = -shortrange1(rj - ri) !!!
!   !!! NOTE THAT phi1(ri - rj) = -phi1(rj - ri) !!!
!
!
!=============================================================================

       SUBROUTINE gammamatrix1(nat,rat,atomtype,basis,period,u, &
                            &  gammamat1x,gammamat1y,gammamat1z)
       IMPLICIT NONE
       INCLUDE 'maxima.inc'
       INTEGER nat,atomtype(*)
       REAL*8 rat(3,*),basis(3,3),tol,u(MAXTP),recbasis(3,3),vol
       REAL*8 gammamat1x(NNDIM,NNDIM),gammamat1y(NNDIM,NNDIM)
       REAL*8 gammamat1z(NNDIM,NNDIM)
       INTEGER i,j
       REAL*8 short_deriv(3),long_deriv(3),r(3),gdrv
       EXTERNAL getalpha
       REAL*8 getalpha,alpha,norm
       LOGICAL period

       
       do i=1,nat
        do j=1,(i-1)
         r(1)=rat(1,j)-rat(1,i)
         r(2)=rat(2,j)-rat(2,i)
         r(3)=rat(3,j)-rat(3,i)
!
!       Lower matrix gammamat1x  contains (phi1+shortrange1)(ri-rj)x
!       Lower matrix gammamat1y  contains (phi1+shortrange1)(ri-rj)y
!       Lower matrix gammamat1z  contains (phi1+shortrange1)(ri-rj)z
!       !!! NOTE THAT gamma1(ri - rj) = -gamma1(rj - ri) !!!
!

       IF (period) THEN

!      get reciprocal lattice vectors and cell volume
       CALL REZVOL(basis,recbasis,vol)

!      choose good convergence parameter alpha
       alpha = getalpha(basis)

!      set tolerance for convergence
       tol = 1d-8
       CALL SHORTRANGE1(r,basis,u(atomtype(j)),u(atomtype(i)), &
                     &  tol,short_deriv)
       CALL PHI1(r,basis,recbasis,alpha,vol,tol,long_deriv)

          gammamat1x(i,j)=short_deriv(1) + long_deriv(1)
          gammamat1y(i,j)=short_deriv(2) + long_deriv(2)
          gammamat1z(i,j)=short_deriv(3) + long_deriv(3)

       ELSE

          norm = sqrt(r(1)**2 + r(2)**2 + r(3)**2)

          call GAM121(norm,u(atomtype(j)),u(atomtype(i)),gdrv)
          gammamat1x(i,j)= gdrv * r(1)
          gammamat1y(i,j)= gdrv * r(2)
          gammamat1z(i,j)= gdrv * r(3)

       END IF

        end do
       end do
       
       END
