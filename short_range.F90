module GWShortRange
  use GWGamma
  implicit none
  private

  public :: shortrange, shortrange1


contains
  !==============================================================================
  !      evaluate short range expression i.e. sumR (gamma - 1/R)
  !      
  !      INPUT:
  !      REAL*8    rh(3)       vector between rmu and rnu
  !      REAL*8    umu         hubbard parameter of orbital mu
  !      REAL*8    unu         hubbard parameter of orbital nu
  !      REAL*8    basis(3,3)  basis of cell(unusual!!!:one line is a basis vector)
  !      REAL*8    tol         convergence tolerance (contribution of last shell)  
  !    
  !      OUTPUT:
  !      REAL*8    value       value of the short range sum
  !==============================================================================

  subroutine SHORTRANGE(rh,basis,umu,unu,tol,value)
    IMPLICIT NONE

    REAL*8 rh(3), umu, unu, basis(3,3), tol, value	
    INTEGER i,j,k,nreal,nmax,nmin 
    REAL*8 rvec(3),R(3),result,lastshell,tmp,norm
    REAL*8 gval

    !
    !      rh = rmu - rnu
    !
    result = 0.0
    nmax = 50
    nmin = 3
    nreal = 0
    lastshell = tol+1d-8
    !      /* sum over R until tolerance is reached */
    DO WHILE ((nreal .le. nmax) .and. ((abs(lastshell) .gt. tol) &
        &  .or. (nreal .le. nmin)  )) 
      lastshell = 0.0
      DO i = -nreal,nreal
        DO j = -nreal,nreal
       	  DO k = -nreal,nreal
            !            /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j))  .or. &
                &  (nreal.eq. abs(k)) ) THEN

              R(1)=i*basis(1,1)+j*basis(2,1)+k*basis(3,1) 
              R(2)=i*basis(1,2)+j*basis(2,2)+k*basis(3,2) 
              R(3)=i*basis(1,3)+j*basis(2,3)+k*basis(3,3) 

              rvec(1) = rh(1) - R(1)  
              rvec(2) = rh(2) - R(2)   
              rvec(3) = rh(3) - R(3)  

              norm   = sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

              !              get value for Gamma
              call GAM12(norm,umu,unu,gval) 

              !              subtract long range part 1/R and multiply by Z(nu)

              IF (norm .lt. 1.0d-6) THEN 
                tmp =  gval
              ELSE
                tmp =  ( gval  - 1/norm )
              ENDIF

              result = result + tmp
              lastshell = lastshell + tmp

            END IF
          END DO
        END DO
      END DO
      nreal = nreal + 1
    END DO

    IF(abs(lastshell) .gt. tol) THEN
      STOP "tolerance in subroutine short not reached."
    END IF
    value = result 

  END subroutine SHORTRANGE




  !=============================================================================
  ! evaluate derivative of short range expression: sumR (d gamma/dR - (-1/R^2))
  !      
  !      INPUT:
  !      REAL*8    rh(3)       vector between rmu and rnu
  !      REAL*8    umu         hubbard parameter of orbital mu
  !      REAL*8    unu         hubbard parameter of orbital nu
  !      REAL*8    basis(3,3)  basis of cell(unusual!!!:one line is a basis vector)
  !      REAL*8    tol         convergence tolerance (contribution of last shell)  
  !    
  !      OUTPUT:
  !      REAL*8    deriv(3)    derivative of the short range sum
  !==============================================================================

  subroutine SHORTRANGE1(rh,basis,umu,unu,tol,deriv)
    REAL*8 rh(3),umu,unu,basis(3,3), tol, deriv(3)
    INTEGER i,j,k,nreal,nmax,nmin
    REAL*8 rvec(3),R(3),lastshell,tmp,norm
    REAL*8 gdrv

    !
    !      rh = rmu - rnu
    !
    deriv(1) = 0.0
    deriv(2) = 0.0
    deriv(3) = 0.0
    nmax = 100
    nmin = 3
    nreal = 0

    lastshell = tol+1.0d-8
    !      /* sum over R until tolerance is reached */
    DO WHILE ((nreal .le. nmax) .and. ( (abs(lastshell) .gt. tol) &
        & .or. (nreal .le. nmin)) ) 
      lastshell = 0.0
      DO i = -nreal,nreal
        DO j = -nreal,nreal
       	  DO k = -nreal,nreal
            !      /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j))  .or. &
                &	(nreal .eq. abs(k)) ) THEN

              R(1)=i*basis(1,1)+j*basis(2,1)+k*basis(3,1) 
              R(2)=i*basis(1,2)+j*basis(2,2)+k*basis(3,2) 
              R(3)=i*basis(1,3)+j*basis(2,3)+k*basis(3,3) 


              rvec(1) = rh(1) - R(1)
              rvec(2) = rh(2) - R(2)
              rvec(3) = rh(3) - R(3)

              norm   = sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

              !       get derivative of gamma 

              call GAM121(norm,umu,unu,gdrv) 

              !       subtract long range -1/R^2/R (see def. of GAM121)
              IF (norm .lt. 1.0d-6) THEN
                tmp = gdrv
              ElSE
                tmp = gdrv - (-1.0/(norm**3)) 
              ENDIF

              deriv(1) = deriv(1) + tmp*rvec(1)
              deriv(2) = deriv(2) + tmp*rvec(2)
              deriv(3) = deriv(3) + tmp*rvec(3)

              lastshell = lastshell + tmp
            END IF
          END DO
        END DO
      END DO
      nreal = nreal + 1
    END DO

    IF(abs(lastshell) .gt. tol) THEN
      STOP "tolerance in subroutine short1 not reached."
    END IF

  end subroutine SHORTRANGE1

end module GWShortRange
