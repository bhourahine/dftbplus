!===============================================================================
!
!    ehfx calculates the HF exchange matrix elements <psi_i|Sigma_x|psi_i>
!    and the DFT exchange correlation matrix elements <psi_i|v_xc|psi_i> 
!    in the basis of Kohn-Sham Orbitals
!
!    The XC part exc is calculated from xc-only SK tables. To account for
!    the density dependence of v_xc, a SCC-correction is applied.
!    The X part is calculated in the gamma-approximation with Hubbard
!    parameters which are derived from Coulomb interaction only. Onsite exchange
!    elements are accounted for additionaly.
!
!   
!    Input: 
!           nn     number of atoms
!           ndim   number of basis functions
!           nang   dimension of GW matrices 
!           ind    indexes AO's belonging to an atom
!           izp    indexes element belonging to an atom
!           rat    atom positions
!           basis  basis of cell  
!           period periodic boundary conditions?
!           c      MO coefficients
!           hxc    XC part of the non-SCC hamiltonian
!           stc    overlap matrix times coefficients
!           occ    occupation numbers
!           dacc   occupation number tolerance
!           qij    transition charges between orbitals i and j
!           gamma  [v] angular momentum resolved, scaled ERI hubbards
!    Output:
!           ehf    <psi_i|Sigma_x|psi_i>
!           exc    <psi_i|v_xc|psi_i> 
!
!===============================================================================
subroutine ehfx (nn,ndim,nang,ind,izp,rat,basis,period, &
               & c,hxc,stc,occ,dacc,qij,gamma,ehf,exc,  &
               & qzero,uhubb,xtab,etab,xm,ntype,lmax) 

  implicit none
  include 'maxima.inc'

  !     Frank-Condon Integrals for carbon
  real*8 g1,f2
  parameter (g1 = 0.267708d0)
  parameter (f2 = 0.173720d0)

  !     Parameter
  integer nn,ndim,nang,ind(*),izp(*)
  real*8  rat(3,*),basis(3,*)
  real*8  c(MDIM,*),hxc(MDIM,*),stc(MDIM,*),occ(*),dacc
  real*8  qij(*),gamma(ADIM,*),ehf(*),exc(*)
  logical period

  !     Local variables
  integer i,j,n,m,k,idm,idn,l,izpk,indl,hocc
  real*8  deltq(NNDIM),gxc(NNDIM,NNDIM)
  real*8  gee(NNDIM,NNDIM),ghu(NNDIM,NNDIM)
  real*8  shift(NNDIM),tmp,tmp1,tmp2
  real*8  cmeri(MAXTYP),cuxc(MAXTYP),intgh(30,30)
  logical indo,ghost

  !     Work arrays
  real*8  wrk1(MDIM,MDIM)

  !     Common block variables
  integer ntype ,lmax(MAXTYP)
  real*8  qzero(MAXTYP)
  real*8  uhubb(MAXTYP) 
  real*8  xtab(MAXTYP,LDIM,LDIM),etab(MAXTYP,LDIM,LDIM)
  real*8  xm(MAXTYP)

  !      indo = .false.
  indo = .true.


  print *,'=====> XC'
  ! ======================== XC part ========================================
  !     Calc <psi_i|v_xc|psi_i>; wrk1 = hxc*c
  call DSYMM('L','U',ndim,ndim,1.d0,hxc,MDIM,c,MDIM,0.d0,wrk1,MDIM)
  do i = 1,ndim
    exc(i) = 0.d0 
    do m = 1,ndim
      exc(i) = exc(i) + c(m,i)*wrk1(m,i)
    end do
  end do

  !     SCC correction to <psi_i|v_xc|psi_i>
  !     dvxc_i = q^i_A * gamma(U_xc)_AB * dq_B; U_xc = U_H - U_ee
  !     Read charges from file
  open(55,file='CHR.DAT')
  do i = 1,5
    read(55,*)
  end do
  do i = 1,nn
    read(55,*) k,deltq(k)
    deltq(k) = deltq(k)-qzero(izp(k))
  end do

  !     Calculate ERI mean over angular momentum
  do i = 1,ntype
    cmeri(i) = 0.d0 
    do l = 1,lmax(i)**2
      do k = l,lmax(i)**2
        cmeri(i) = cmeri(i) +  etab(i,l,k)
      end do
    end do
    cmeri(i) = 2.d0*cmeri(i)/(lmax(i)**4+lmax(i)**2)
    cuxc(i)=-(uhubb(i)-cmeri(i))
  end do
  !      call gammamatrix(nn,rat,izp,basis,period,cuxc,gxc)
  call gammamatrix(nn,rat,izp,basis,period,uhubb,ghu)
  call gammamatrix(nn,rat,izp,basis,period,cmeri,gee)
  do i = 1,nn
    do j = 1,nn
      gxc(i,j) = ghu(i,j) - gee(i,j)
    end do
  end do

  print *,'=====> SCC-XC'
  !     shift = gxc*deltq
  call DSYMV('L',nn,1.d0,gxc,NNDIM,deltq,1,0.d0,shift,1)
  do i = 1,ndim
    call transq(nn,stc,c,qij,i,i,lmax,izp)
    indl = 0
    do k = 1,nn
      izpk = izp(k)
      tmp = 0.d0            
      do l = 1,lmax(izpk)
        indl = indl + 1
        tmp = tmp +  qij(indl)
      end do
      exc(i) = exc(i) +  tmp*shift(k)
    end do
  end do
  print *,'=====> HF'
  ! 
  ! ======================== X part ======================================== 
  !     Calculate exchange energy per orbital CNDO

  !     Find highest occupied orbital
  do k = 1,ndim
    if(abs(occ(k)).lt.dacc) exit
  end do
  hocc = k-1

  do i = 1,ndim
    ehf(i) = 0.0d0
    do k = 1,hocc
      call transq(nn,stc,c,qij,i,k,lmax,izp)
      do n = 1,nang
        do m = 1,nang
          ehf(i) = ehf(i) - &
              & 0.5d0*occ(k)*qij(m)*gamma(m,n)*qij(n)
        end do
      end do
    end do
  end do


  print *,'=====> INDO'
  !     Add one-center exchange terms

  if(indo) then
    !       wrk1 = sum_l^nhel c(m,l)*(n,l)
    do m = 1,ndim
      do n = 1,ndim
        wrk1(n,m) = 0.d0
        do l = 1,hocc
          wrk1(n,m) =  wrk1(n,m) + 0.5d0*occ(l)*c(m,l)*c(n,l)
        end do
      end do
    end do
    do i = 1,ndim
      do k = 1,nn
        do m = 1,ind(k+1)-ind(k)
          idm = ind(k)+m
          do n = m+1,ind(k+1)-ind(k)
            idn = ind(k)+n
            tmp1 = wrk1(idn,idn)*(c(idm,i)**2) &
                & + wrk1(idm,idm)*(c(idn,i)**2)&
                & + 2*wrk1(idm,idn)*c(idm,i)*c(idn,i)
            ehf(i) = ehf(i) - tmp1*xtab(izp(k),m,n)  
          end do
        end do
      end do
    end do

    !       Control if there are ghost atoms
    ghost = .false.
    do i=1,nn
      if ( abs(xm(izp(i))) .lt. 0.1) then
        ghost = .true.
        call ghostint(intgh,lmax)
      end if
    end do

    !       Important, this part contains the code to manage the
    !       exception represented by ghost atoms. It is used to enlarge the basis of metallic
    !       atoms: it is mandatory that every ghost atom in the input geometry file must be
    !       before the respective metallic atom.
    if(ghost) then
      do i=1,ndim
        do k=1,nn
          if ( abs(xm(izp(k))) .lt. 0.1 ) then
            do m = 1, ind(k+1) - ind(k)
              idm = ind(k)+m
              do n = 1,ind(k+2) - ind(k+1)
                idn = ind(k+1) + n
                !      Because the ghost is just before the metallic atom the loop must run 
                !      between ind(k) and ind(k+1) for the ghost and ind(k+1) and ind(k+2)
                !      for the metallic atom.
                tmp = wrk1(idn,idn)*(c(idm,i)**2) &
                  & + wrk1(idm,idm)*(c(idn,i)**2) &
                  & + 2*wrk1(idm,idn)*c(idm,i)*c(idn,i)
                ehf(i) = ehf(i) - tmp*intgh(m,n)
              end do
            end do
          end if
        end do
      end do
    end if
  end if
  return
end subroutine ehfx

