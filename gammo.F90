module GWGammo
  use SimpleAlgebra
  use GWLongRange
  use GWShortRange
  use GWGamma
  implicit none
  private

  public :: gammo

contains


  !=============================================================================
  !   Build lower triangular matrix containing Ewald potentials + short
  !   range terms
  !
  !   INPUT Parameter:
  !   INTEGER nat          number of atoms
  !   REAL*8 rat(3,*)      position of atoms
  !   REAL*8 u(MAXTYP,*)   hubbard parameters, l resolved
  !   REAL*8 basis(3,3)    basis of cell
  !   LOGICAL period       .true. if periodic boundary conditions
  !
  !   OUTPUT:
  !   REAL*8 gammamat(ADIM,*) matrix containing the values of the ewlad potential
  !                       in the upper triangular part
  !      !!! NOTE THAT phi(ri - rj) = phi(rj - ri) !!!
  !      !!! NOTE THAT shortrange(ri - rj) = shortrange(rj - ri) !!!
  !
  !=============================================================================

  SUBROUTINE gammo(nat,rat,atomtype,basis,period,u,lind,lmax, &
      & gammamat,nType,nAng)
    INTEGER nat
    REAL*8 rat(3,*),basis(3,3),tol,u(nType,*),gammamat(nAng,nAng)
    INTEGER i,j,atomtype(*),lind(*),lmax(*)
    INTEGER izpi,izpj,lofi,lofj
    REAL*8 tmpvalue,phivalue,r(3)
    REAL*8 recbasis(3,3), vol
    REAL*8 alpha
    REAL*8 gval,norm
    LOGICAL period

    integer, intent(in) :: nType
    integer, intent(in) :: nAng

    do i=1,nat
      izpi=atomtype(i)
      do j=1,i
        izpj=atomtype(j)
        r(1)=rat(1,i)-rat(1,j)
        r(2)=rat(2,i)-rat(2,j)
        r(3)=rat(3,i)-rat(3,j)

        !       !!! NOTE THAT phi(ri - rj) = phi(rj - ri) !!!
        !       !!! NOTE THAT shortrange(ri - rj) = shortrange(rj - ri) !!!


        IF(period) THEN
          vol = determinant33(basis)
          call invert33(recbasis, basis, vol)
          recbasis = reshape(recbasis, (/3, 3/), order = (/2, 1/))
          !      choose good convergence parameter alpha
          alpha = getalpha(basis)
          !      set tolerance for convergence
          tol = 1.0d-8

          CALL phi(r,basis,recbasis,alpha,vol,tol,phivalue)
          do lofi=1,lmax(izpi)
            do lofj=1,lmax(izpj)
              CALL shortrange(r,basis,u(izpi,lofi),u(izpj,lofj), &
                  &  tol,tmpvalue)
              gammamat(lind(i)+lofi,lind(j)+lofj)= &
                  &  tmpvalue + phivalue    
            enddo
          enddo

        ELSE

          norm   = sqrt(r(1)**2+r(2)**2+r(3)**2)
          do lofi=1,lmax(izpi)
            do lofj=1,lmax(izpj)
              !      get value for Gamma
              CALL GAM12(norm,u(izpi,lofi),u(izpj,lofj),gval)
              gammamat(lind(i)+lofi,lind(j)+lofj)=gval
            enddo
          enddo


        END IF

      END DO
    END DO

  END SUBROUTINE gammo

end module GWGammo
