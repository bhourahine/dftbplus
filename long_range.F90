!  =====================================================================
!  evaluation of the potential phi 
!
!   phi(r) = 4*pi/Omega ( Summe{G neq 0} e^{-G^2/{4 alpha^2}}/{G^2} cos(G r)
!            +Summe{R, R neq r} (1-erf(alpha*|R-r|))/|R-r|
!            -pi/(Omega*alpha^2)
!
!
!   INPUT Parameter:
!   REAL*8 r(3)           position of evaluation of potential
!   REAL*8 basis(3,3)     basis of cell
!   REAL*8 recbasis(3,3)      basis of reciprocal cell
!   REAL*8 alpha          convergence parameter
!   REAL*8 vol            cell volume
!   REAL*8 tol            tolerance for convergence of "last shell"
!                       (can often be used as a criterion for global
!                        convergence) 
!   OUTPUT:
!   REAL*8 potential      value of Ewald potential
!  
!  ======================================================================

        subroutine phi(r,basis,recbasis,alpha,vol,tol,potential)
        IMPLICIT NONE
        external terfc
        REAL*8 terfc
        REAL*8 r(3), basis(3,3), recbasis(3,3), alpha, vol, tol
        REAL*8 potential
        REAL*8 reciprocal,rspace,cterm
        REAL*8 G(3),rh(3),help,norm,lastshell
        REAL*8 MPI
        INTEGER nrezi, nreal, nmax, nmin
        INTEGER i,j,k
        MPI = 3.14159265358979323846 
        nmax = 20
        nmin = 2

!       evaluate reciprocal space term ( sum over G <> 0) ...  
!       /* sum over G until tolerance is reached */
        nrezi = 1
        lastshell = tol+1d-8  
        reciprocal = 0.0
        DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or. &
        & (abs(lastshell) .gt.  tol)))
          lastshell = 0.0  
          DO i=-nrezi,nrezi
           DO j=-nrezi,nrezi
            DO k=-nrezi,nrezi
!             /*only G belonging to outer shells are new ones */
              IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or. &
            & (nrezi. eq. abs(k)) ) THEN
                  G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
                  G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
                  G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)
 
                  help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
                  help = exp(-help/(4.0*alpha*alpha))/help
                  help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help
 
                  reciprocal = reciprocal + help
                  lastshell = lastshell + help/vol       
              END IF
             END DO
            END DO
          END DO
         nrezi = nrezi + 1
        END DO
     
!       stop if tolerance not reached
        IF ( abs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi not reached in reciprocal space"     
        END IF

        reciprocal=(4.0*MPI*reciprocal)/vol


!       evaluate  real space term (sum over R)   
!       /* sum over R until tolerance is reached */
	rspace = 0.0
        nreal = 0
         lastshell = tol+1d-8  
        DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin) & 
                & .or. (abs(lastshell) .gt.  tol)))
         lastshell = 0.0  
         DO i=-nreal,nreal
          DO j=-nreal,nreal
           DO k=-nreal,nreal
!            /*only R belonging to outer shells are new ones */
             IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or. &
             & (nreal .eq. abs(k)) ) THEN
                rh(1)=r(1)-(i*basis(1,1)+j*basis(2,1)+k*basis(3,1))
                rh(2)=r(2)-(i*basis(1,2)+j*basis(2,2)+k*basis(3,2))
                rh(3)=r(3)-(i*basis(1,3)+j*basis(2,3)+k*basis(3,3))
                norm=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3))
                IF (norm .gt. 1.0d-20) THEN 
!               erfc=1-erf   
                  help   = terfc(alpha*norm)/norm
                  rspace = rspace + help
                  lastshell = lastshell + help
                ELSE 
                   lastshell = tol+1d-8  
              END IF
             END IF
            END DO
           END DO
          END DO
         nreal = nreal + 1
        END DO

!       stop if tolerance not reached
        IF ( abs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi not reached in real space"     
        END IF


!       evaluate constant term pi/(Omega*alpha^2) 
        cterm = -MPI/(vol*alpha*alpha)

!       if r = 0 there is another constant to be added   	
        IF ((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) .lt. 1.0d-20) THEN
            cterm = cterm -2.0*alpha/sqrt(MPI)
        END IF

!       set the value of the potential
        potential = reciprocal + rspace + cterm

        END 


 
!  =====================================================================
!  evaluation of the derivative of the potential phi 
!
!   INPUT Parameter:
!   REAL*8 r(3)           position of evaluation of potential
!   REAL*8 basis(3,3)     basis of cell
!   REAL*8 recbasis(3,3)      basis of reciprocal cell
!   REAL*8 alpha          convergence parameter
!   REAL*8 vol            cell volume
!   REAL*8 tol            tolerance for convergence of "last shell"
!                       (can often be used as a criterion for global
!                        convergence) 
!   OUTPUT:
!   REAL*8 deriv(3)       derivative of ewlad potential
!  
!  ======================================================================

 
	subroutine phi1(r,basis,recbasis, alpha,vol,tol,deriv)

        IMPLICIT NONE
        external terfc
        REAL*8 terfc 
        REAL*8 r(3), basis(3,3), recbasis(3,3), alpha, vol, deriv(3)
        REAL*8 reciprocal(3),rspace(3),MPI 
        REAL*8 G(3),rh(3),norm,help,tol,lastshell 
        INTEGER i,j,k, nrezi, nreal, nmax, nmin 
        nmax = 20
        nmin = 2
        MPI =  3.14159265358979323846

        IF((r(1)*r(1) + r(2)*r(2) + r(3)*r(3)) .lt. 1.0d-20) THEN
          deriv(1) = 0.0
          deriv(2) = 0.0
          deriv(3) = 0.0
          return
        END IF
 
!       /* evaluate reciprocal space term (sum over G <> 0) ...  */
        nrezi = 1
        lastshell = tol+1d-8
       	reciprocal(1) = 0.0  
        reciprocal(2) = 0.0  
        reciprocal(3) = 0.0 
        DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or. &
        & (abs(lastshell) .gt.  tol)))
          lastshell = 0.0
        DO i=-nrezi,nrezi
         DO j=-nrezi,nrezi
          DO k=-nrezi,nrezi
!             /*only G belonging to outer shells are new ones */
              IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or. &
              & (nrezi. eq. abs(k)) ) THEN
 
                  G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1) 
                  G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2) 
   	          G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3) 
 
   	          help=G(1)*G(1)+G(2)*G(2)+G(3)*G(3) 
 
   	          help=exp(-help/(4.0*alpha*alpha))/help 
 
   	          help=-sin(G(1)*r(1)+G(2)*r(2)+G(3)*r(3))*help 
 
   	          reciprocal(1)=help*G(1) + reciprocal(1)
   	          reciprocal(2)=help*G(2) + reciprocal(2)
   	          reciprocal(3)=help*G(3) + reciprocal(3)

                  lastshell = lastshell + help/vol
              ENDIF
          END DO
         END DO
        END DO
        nrezi = nrezi + 1
        END DO
 
!       stop if tolerance not reached
        IF ( abs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi1 not reached in reciprocal space"
        END IF

        reciprocal(1)=(4.0*MPI*reciprocal(1))/vol 
        reciprocal(2)=(4.0*MPI*reciprocal(2))/vol 
        reciprocal(3)=(4.0*MPI*reciprocal(3))/vol 
 
 
!       /* evaluate  real space term (sum over R) */
!       /* sum over R until tolerance is reached */
        rspace(1) = 0.0  
        rspace(2) = 0.0  
        rspace(3) = 0.0 
        nreal = 0
        lastshell = tol+1d-8
        DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin) &
                 & .or. (abs(lastshell) .gt.  tol)))
        lastshell = 0.0
        DO i=-nreal,nreal
         DO j=-nreal,nreal
          DO k=-nreal,nreal
!            /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or. &
           & (nreal. eq. abs(k)) ) THEN
            rh(1)=r(1)-(i*basis(1,1)+j*basis(2,1)+k*basis(3,1)) 
            rh(2)=r(2)-(i*basis(1,2)+j*basis(2,2)+k*basis(3,2)) 
            rh(3)=r(3)-(i*basis(1,3)+j*basis(2,3)+k*basis(3,3)) 
 
            norm=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3)) 
 
            help = (-2/sqrt(MPI)*exp(-alpha*alpha*norm*norm)* &
                 & alpha*norm - terfc(alpha*norm))/(norm*norm*norm) 
 
            rspace(1) = rh(1)*help + rspace(1)  
            rspace(2) = rh(2)*help + rspace(2)
            rspace(3) = rh(3)*help + rspace(3)
   
            lastshell = lastshell + help
            END IF
          END DO
         END DO
        END DO
       nreal = nreal + 1
      END DO

!       stop if tolerance not reached
        IF ( abs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi1 not reached in real space"
        END IF


!       /* add real and reciprocal parts */
        deriv(1) = rspace(1)  + reciprocal(1) 
        deriv(2) = rspace(2)  + reciprocal(2) 
        deriv(3) = rspace(3)  + reciprocal(3) 

        END 



!==============================================================================
!
!       evaluate the cross product of A x B
!
!==============================================================================

	subroutine CROSS( A, B, C) 
	IMPLICIT NONE
	REAL*8 A(3), B(3), C(3)

      	  C(1)=A(2)*B(3)-A(3)*B(2)
      	  C(2)=A(3)*B(1)-A(1)*B(3)
     	  C(3)=A(1)*B(2)-A(2)*B(1)
	END 


!       get reciprocal lattice vectors and volume of unit cell       

	subroutine REZVOL(basis,recbasis,vol) 
	IMPLICIT NONE
	external cross
	REAL*8  MPI
	REAL*8  basis(3,3), recbasis(3,3), vol
	REAL*8  hv1(3), hv2(3), hv3(3), hv4(3), fac
       	INTEGER i
	MPI = 3.14159265358979323846 
		

       	  DO i=1,3, 1
            hv1(i)=basis(1,i)
            hv2(i)=basis(2,i)
            hv3(i)=basis(3,i)
	  END DO

	  call CROSS(hv2,hv3,hv4)
	  vol = abs(basis(1,1)*hv4(1)+basis(1,2)*hv4(2)+basis(1,3)*hv4(3))
	  fac = 2.0*MPI/vol

	  recbasis(1,1)=hv4(1)*fac
	  recbasis(1,2)=hv4(2)*fac
	  recbasis(1,3)=hv4(3)*fac

	  call CROSS(hv3,hv1,hv4)
	  recbasis(2,1)=hv4(1)*fac 
	  recbasis(2,2)=hv4(2)*fac 
	  recbasis(2,3)=hv4(3)*fac

	  call CROSS(hv1,hv2,hv4)
       	  recbasis(3,1)=hv4(1)*fac
	  recbasis(3,2)=hv4(2)*fac
	  recbasis(3,3)=hv4(3)*fac
		
	END 



!==============================================================================
!
!     Returns the (tabulated) complementary error function erfc(x) 
!     with fractional error everywhere less than 1.2 x 10^-7
!
!==============================================================================

	FUNCTION terfc(x)
	REAL*8 terfc,x
	REAL*8 t,z

	z = abs(x)
	t = 1.0/(1.0+0.5*z)
	terfc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+ &
       & t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+ &
       & t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
	if (x .lt. 0.0) terfc=2.0-terfc
	
	RETURN
	END


!==============================================================================
!      get optimal alpha for Ewald potential phi
!
!      INPUT:
!      REAL*8 basis(3,3)     basis of lattice                     
!
!      RETURNS:
!      REAL*8                optimal alpha
!==============================================================================

       function getalpha(basis)
       IMPLICIT NONE
        REAL*8 basis(3,3)
        REAL*8 getalpha
        REAL*8 alpha, alphal, alphar
        INTEGER nopt
        REAL*8 recbasis(3,3), vol, tol
        REAL*8 G, R, help1, help2, help3
        EXTERNAL diffrecreal
        REAL*8 diffrecreal
       
        tol  = 1d-5

!       get reciprocal lattice vectors and cell volume

        CALL REZVOL(basis,recbasis,vol)

!       get sqnorm of smallest vector in reciprocal space 
        help1 = recbasis(1,1)**2+recbasis(1,2)**2+recbasis(1,3)**2
        help2 = recbasis(2,1)**2+recbasis(2,2)**2+recbasis(2,3)**2
        help3 = recbasis(3,1)**2+recbasis(3,2)**2+recbasis(3,3)**2
        G    = sqrt(min(help1,help2,help3))
        
!       get norm of smallest vector in real space 
        help1 = basis(1,1)**2 + basis(1,2)**2 + basis(1,3)**2
        help2 = basis(2,1)**2 + basis(2,2)**2 + basis(2,3)**2
        help3 = basis(3,1)**2 + basis(3,2)**2 + basis(3,3)**2
        R    = sqrt(min(help1,help2,help3))


!       optimise alpha

!       if (reciprocalspace decline - realspace decline) < 0 convergence
!       in real space too slow: increase alpha
!       set starting alphal
        alpha = 1d-5
         DO WHILE( diffrecreal(alpha,G,R,vol) .lt. tol )
          alphal = alpha
          alpha = alpha*2.0   
         END DO

!       if (reciprocalspace decline - realspace decline) > 0 convergence
!       in reciprocal space too slow: decrease alpha
!       set starting alphar
        alpha = 1e+5
         DO WHILE( diffrecreal(alpha,G,R,vol) .gt. tol )
          alphar = alpha
          alpha = alpha / 2.0
         END DO

!       now find best value by refining the interval [alphal,alphar]
         alpha = (alphal + alphar)/2.0
         nopt  = 0
        DO WHILE ( abs(diffrecreal(alpha,G,R,vol)) .gt. tol .and. &
                & nopt .le. 20)
         IF ( diffrecreal(alpha,G,R,vol) .lt. tol ) THEN
          alphal = alpha
         END IF

         IF ( diffrecreal(alpha,G,R,vol) .gt. tol ) THEN
          alphar = alpha
         END IF
         alpha = (alphal + alphar)/2.0
         nopt  = nopt + 1
        END DO

        IF (nopt .gt. 20 ) THEN 
        alpha=exp(-0.310104*log(vol)+0.786382)/2.0
         PRINT*, "WARNING: NO OPTIMISED ALPHA FOUND: "
         PRINT*, "STANDARD ALPHA USED. ALPHA SET TO", alpha 
        END IF
          
        getalpha = alpha

        END


!==============================================================================
!      
!      This subroutine returns the norm of the largest reciprocal space
!      vector and the norm of the largest real space vector needed for 
!      a converged Ewald summation
!
!      INPUT:
!      REAL*8  alpha         chosen convergence parameter
!      REAL*8  tol           disired tolerance of last contribution
!      REAL*8  vol           cell volume
!
!      OUTPUT:
!      REAL*8  Gmax          norm of largest vector in reciprocal space
!      REAL*8  Rmax          norm of largest vector in real space 
!==============================================================================


       subroutine getGRmax(alpha,tol,vol,Gmax,Rmax)
       IMPLICIT NONE
        REAL*8 alpha, Gmax, Rmax, tol, vol
        INTEGER nopt
        REAL*8 G, R, Gl, Gr, Rl, Rr
        EXTERNAL Gspace 
        REAL*8 Gspace 
        EXTERNAL Rspace 
        REAL*8 Rspace 
        EXTERNAL REZVOL
        


!       set starting Gl and Rl
        G = 1d-5
        R = 1d-5
         DO WHILE( Gspace(G,alpha,vol) .gt. tol + 1d-10)
          Gl = G
          G    = G*2.0   
         END DO

         DO WHILE( Rspace(R,alpha) .gt. tol + 1d-10)
          Rl = R
          R    = R*2.0   
         END DO

!       set starting Gr and Rr
        G = 1e+5
        R = 1e+5
         DO WHILE( Gspace(G,alpha,vol) .lt. tol - 1d-10)
          Gr = G
          G    = G/2.0   
         END DO

         DO WHILE( Rspace(R,alpha) .lt. tol - 1d-10)
          Rr = R
          R    = R/2.0   
         END DO



!       now find best value of G by refining the interval [Gl,Gr]
          G = (Gl + Gr)/2.0
         nopt  = 0
        DO WHILE ( Gspace(G,alpha,vol) .gt. tol .and. nopt .le. 20)
         IF ( Gspace(G,alpha,vol) .lt. (tol - 1d-10) ) THEN
          Gr = G
         END IF

         IF ( Gspace(alpha,G,vol) .gt. (tol + 1d-10) ) THEN
          Gl = G
         END IF

         G = (Gl + Gr)/2.0
         nopt  = nopt + 1
        END DO

        IF (nopt .ge. 20) THEN
          STOP "Gmax COULD NOT BE DETERMINED IN getGRmax"
        END IF

!       now find best value of R by refining the interval [Rl,Rr]
          R = (Rl + Rr)/2.0
         nopt  = 0
        DO WHILE ( Rspace(R,alpha) .gt. tol .and. nopt .le. 20)

         IF ( Rspace(R,alpha) .lt. (tol - 1d-10) ) THEN
          Rr = R
         END IF

         IF ( Rspace(alpha,R) .gt. (tol + 1d-10) ) THEN
          Rl = R
         END IF

         R = (Rl + Rr)/2.0
         nopt  = nopt + 1
        END DO

        IF (nopt .ge. 20) THEN
          STOP "Rmax COULD NOT BE DETERMINED IN getGRmax"
        END IF

        Gmax = G
        Rmax = R
        PRINT*,"GMAX:", GMAX, "REC:", Gspace(Gmax,alpha,vol)
        PRINT*,"RMAX:", RMAX, "REAL*8:", Rspace(Rmax,alpha)

        END



!==============================================================================
!      get differnce between decline in reciprocal and real space
!      this function is only used by function getalpha
!
!      INPUT:
!      REAL*8  alpha         convergence parameter
!      REAL*8  G             square of norm of smallest G
!      REAL*8  R             norm of smallest R
!      REAL*8  vol           cell volume
!
!      RETURNS:
!      REAL*8                difference between decline in reciprocal 
!                          space (rec(2G)-rec(3G)) and real space (real(2R) 
!                           - real(3R))
!==============================================================================

        function diffrecreal(alpha,G,R,vol)
        IMPLICIT NONE
        REAL*8 alpha, G, R, vol
        REAL*8 diffrecreal
        REAL*8 diffrec, diffreal
        EXTERNAL Gspace
        REAL*8 Gspace
        EXTERNAL Rspace
        REAL*8 Rspace
        REAL*8 MPI
        MPI = 3.14159265358979323846 

!       make differences between decline at 2G and 3G / 2R and 3R
        diffrec = Gspace(2.0*G,alpha,vol) - Gspace(3.0*G,alpha,vol)
        diffreal= Rspace(2.0*R,alpha) - Rspace(3.0*R,alpha)

!       return difference between reciprocal and realspace decline
        diffrecreal        = diffrec - diffreal

        END


!==============================================================================
!       returns the "r independent" G space part of the Ewald sum
!      
!       INPUT:
!       REAL*8 G       norm of G
!       REAL*8 alpha   chosen convergence parameter
!       REAL*8 vol     cell volume
!
!       RETURNS:     
!       REAL*8         "r independent" G space part of the Ewald sum
!============================================================================== 

        function Gspace(G,alpha,vol)
 
        IMPLICIT NONE
        REAL*8 G, alpha, vol
        REAL*8 Gspace
        REAL*8 MPI
        MPI = 3.14159265358979323846 

!       evaluate reciprocal space term at G
        Gspace = exp(-G**2/(4.0*alpha*alpha))/(G**2)
        Gspace = (4.0*MPI*Gspace)/vol

        END


!==============================================================================
!       returns the R space part of the Ewald sum
!      
!       INPUT:
!       REAL*8 R       norm of R
!       REAL*8 alpha   chosen convergence parameter
!
!       RETURNS:     
!       REAL*8         R space part of the Ewald sum
!============================================================================== 

        function Rspace(R,alpha)

        IMPLICIT NONE
        REAL*8 R, alpha
        REAL*8 Rspace
        EXTERNAL terfc
        REAL*8 terfc

!       evaluate real space term at R
        Rspace  = terfc(alpha*R)/R

        END
