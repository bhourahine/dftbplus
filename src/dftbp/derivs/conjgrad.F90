!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Conjugate gradient solvers for problems of the forms Ab = x
module dftbp_derivs_conjgrad
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_sparseblas, only : symv, symm
  use dftbp_io_message, only : warning
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: conjgrad, precondCG ! conjgrad_contract, bicgstab

  !> Calculates b in A b = x using halting criteria |A b - x| < tol
  interface conjgrad
    module procedure dense
    module procedure sparse
    module procedure sparse_contract
  end interface conjgrad

  !> Calculates b in A b = x using a pre-conditioning matrix
  interface precondCG
    module procedure densePreCond
  end interface precondCG

!  !> Calculates b in A b = x for complex A
!  interface bicgstab
!    module procedure bicgstab_cmplx
!  end interface bicgstab

  character(lc) :: lcTmp

contains

  !> Calculates b in Ab = x using halting criteria |A b - x| < tol
  !> note based on wikipedia pseudo code, no pre-conditioning
  subroutine dense(A, b, x, tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> solution vector
    real(dp), intent(in) :: b(:)

    !> on entry intial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), p(size(b)), Ap(size(b))
    real(dp) :: rsold, rsnew, alpha, tol_
    integer :: ii

    @:ASSERT(size(A,dim = 1) ==  size(A,dim = 2))
    @:ASSERT(size(b) ==  size(x))
    @:ASSERT(size(A,dim = 1) ==  size(x))

    tol_ = setR2Tol(tol)

    r(:) = b - matmul(A, x)
    p(:) = r
    rsold = dot_product(r, r)
    do ii = 1, size(A, dim=1)
      Ap(:) = matmul(A, p)
      alpha = rsold / dot_product(p, Ap)
      x(:) = x + alpha * p
      r(:) = r - alpha * Ap
      rsnew = dot_product(r, r)
      if (rsnew < tol_) then
        exit
      end if
      p(:) = r + (rsnew / rsold) * p
      rsold = rsnew
    end do

  end subroutine dense


  !> Calculates b in Ab = x using halting criteria |A b - x| < tol
  !> note based on wikipedia pseudo code, pre-conditioning by M
  subroutine densePreCond(A, Minv, b, x, tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> inverse of pre-conditioner matrix M
    real(dp), intent(in) :: Minv(:,:)

    !> solution vector
    real(dp), intent(in) :: b(:)

    !> on entry intial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: rNew(size(b)), zNew(size(b)), p(size(b)), Ap(size(b))
    real(dp) :: rOld(size(b)), zOld(size(b))
    real(dp) :: rsOld, rsNew, alpha, beta, tol_
    integer :: ii

    @:ASSERT(size(A,dim = 1) ==  size(A,dim = 2))
    @:ASSERT(size(b) ==  size(x))
    @:ASSERT(size(A,dim = 1) ==  size(x))

    tol_ = setR2Tol(tol)

    rOld(:) = b - matmul(A, x)
    zOld(:) = matmul(Minv, rOld)
    p(:) = zOld
    rsOld = dot_product(rOld, rOld)
    do ii = 1, size(A, dim=1)
      Ap(:) = matmul(A, p)
      alpha = rsOld / dot_product(p, Ap)
      x(:) = x + alpha * p
      rNew(:) = rOld - alpha * Ap
      rsNew = dot_product(rNew, rNew)
      if (rsNew < tol_) then
        exit
      end if
      zNew(:) = matmul(Minv, rNew)
      ! Polak-Ribière
      beta = dot_product(rNew, zNew - zOld) / dot_product(rOld, zOld)
      p(:) = zNew + beta * p
      rsOld = rsNew
      rOld(:) = rNew
      zOld(:) = zNew
    end do

  end subroutine densePreCond


  !> Version of the dense subroutine using the DFTB+ sparse matrix stored A and calculating b in
  !> Ab = x using halting criteria |A b - x| < tol, no pre-conditioning
  subroutine sparse(A, b, x, iNeighbor, nNeighbor, img2CentCell, iSparseStart, iAtomStart, orb, tol)

    !> Symmetric positive definite matrix in sparse storage
    real(dp), intent(in) :: A(:)

    !> Solution vector
    real(dp), intent(in) :: b(:)

    !> On entry intial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), p(size(b)), Ap(size(b)), rsold, rsnew, alpha, tol_
    integer :: ii

    @:ASSERT(size(b) ==  size(x))

    tol_ = setR2Tol(tol)

    ! replaces r = b - matmul(A, x)
    r(:) = b
    call symv(r, A, x, iNeighbor, nNeighbor, img2CentCell, iSparseStart, iAtomStart, orb,&
        & alpha=-1.0_dp, beta=1.0_dp)
    p(:) = r
    rsold = dot_product(r,r)
    do ii = 1, size(A,dim = 1) ! max iterations to coverge
      ! replaces Ap = matmul(A,p);
      call symv(Ap, A, p, iNeighbor, nNeighbor, img2CentCell, iSparseStart, iAtomStart, orb)
      alpha = rsold / dot_product(p,Ap)
      x(:) = x + alpha * p
      r(:) = r - alpha * Ap
      rsnew = dot_product(r, r)
      if (rsnew < tol_) then
        exit
      end if
      p(:) = r + (rsnew / rsold) * p
      rsold = rsnew
    end do

  end subroutine sparse


  !> Sparse version multiple RHS CG (no pre-conditioning) for solving A B = X, or A (B - eS) = X if
  !> overlap matrix provided. Reduces iterating vectors as they converge
  subroutine sparse_contract(A, B, X, iNeighbor, nNeighbor, img2CentCell, iSparseStart, iAtomStart,&
      & orb, ei, over, tol, nCullIter)

    !> Symmetric positive definite matrix in sparse storage
    real(dp), intent(in) :: A(:)

    !> Solution vectors
    real(dp), intent(in) :: B(:,:)

    !> On entry intial guesses to solutions, on exit the solutions
    real(dp), intent(inout) :: X(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Eigenvalues
    real(dp), intent(in), optional :: ei(:)

    !> Overlap in sparse format
    real(dp), intent(in), optional :: over(:)

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    !> Period at which to remove converged vectors from iteration
    integer, intent(in), optional :: nCullIter

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xTmp(:,:)
    logical, allocatable :: isNotConverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u, nCullIter_
    integer, allocatable :: indx(:)

    @:ASSERT(present(ei).eqv.present(over))

    tol_ = setR2Tol(tol)
    nCullIter_ = setIVal(40, nCullIter)
    @:ASSERT(nCullIter_ > 0)

    n = size(x,dim = 1)
    s = size(x,dim = 2)
    @:ASSERT(all(shape(b) == [n,s]))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(isNotConverged(s))
    allocate(indx(s))
    if (present(over)) then
      allocate(xTmp(n,s))
    end if

    t = s
    do ii = 1, t
      indx(ii) = ii
    end do

    isNotConverged(:) = .true.

    r(:,:) = B
    call symm(r, 'L', A, x, iNeighbor, nNeighbor, img2CentCell, iSparseStart, iAtomStart, orb,&
        & alpha=-1.0_dp, beta=1.0_dp)
    if (present(over)) then
      call symm(xTmp, 'L', over, x, iNeighbor, nNeighbor, img2CentCell, iSparseStart, iAtomStart,&
          & orb)
      do jj = 1, t
        r(:,jj) = r(:,jj) + ei(jj) * xTmp(:,jj)
      end do
    end if
    p(:,:) = r
    rsOld(:) = sum(r**2, dim=1)

    do ii = 1, size(B, dim=1)

      call symm(Ap(:,:t), 'L', A, p(:,:t), iNeighbor, nNeighbor, img2CentCell, iSparseStart,&
          & iAtomStart, orb)
      if (present(over)) then
        call symm(xTmp(:,:t), 'L', over, p(:,:t), iNeighbor, nNeighbor, img2CentCell, iSparseStart,&
            & iAtomStart, orb)
        do jj = 1, t
          Ap(:,jj) = Ap(:,jj) - ei(indx(jj)) * xTmp(:,jj)
        end do
      end if
      alpha(:t) = rsold(:t) / sum(p(:,:t) * Ap(:,:t),dim=1)

      do jj = 1, t
        if (isNotConverged(jj)) then
          x(:,indx(jj)) = x(:,indx(jj)) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsNew(:t) = sum(r(:,:t)**2,dim = 1)

      do jj = 1, t
        if (isNotConverged(jj) .and. rsNew(jj)<tol_) then
          isNotConverged(jj) = .false.
        end if
      end do

      if (t > 1 .and. mod(ii, nCullIter_) ==  0) then
        ! contract list of vectors to solve
        u = 1
        do jj = 1, t
          if (isNotConverged(jj)) then
            if (u < jj) then
              indx(u) = indx(jj)
              isNotConverged(u) = isNotConverged(jj)
              Ap(:,u) = Ap(:,jj)
              p(:,u) = p(:,jj)
              r(:,u) = r(:,jj)
              rsnew(u) = rsnew(jj)
              rsold(u) = rsold(jj)
              alpha(u) = alpha(jj)
            end if
            u = u +1
          end if
        end do
        t = u - 1
      end if
      if (all(.not.isNotConverged(:t)) .or. t < 1) then
        exit
      end if

      do jj = 1, t
        p(:,jj) = r(:,jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:t) = rsnew(:t)

    end do

    if (any(isNotConverged)) then
      write(lcTmp,"('Not all vectors converged in CG ',I0, ' failed')") count(isNotConverged)
      call warning(trim(lcTmp))
    end if

  end subroutine sparse_contract


!  !> Modified version of BICGSTAB algorithm to solve Ax = b, where the interation
!  !! restarts itself if alpha goes to 0 (so if it were about to fail on the
!  !! next iteration).
!  !! \param x on entry intial guess to solution, on exit the solution
!  !! \param a matrix
!  !! \param b right hand side
!  !! \param tol tolerance for | Ax - b |
!  subroutine bicgstab_cmplx(a,b,x,tol)
!    complex(dp), intent(in) :: A(:,:)
!    complex(dp), intent(in) :: b(:)
!    complex(dp), intent(inout) :: x(:)
!    real(dp), intent(in) :: tol
!
!    complex(dp), allocatable :: r(:,:), rh(:,:), nu(:,:), p(:,:), xtmp(:,:)
!    complex(dp), allocatable :: rho(:), omega(:)
!    complex(dp) :: s(size(b)), t(size(b)), tmp(size(b))
!    complex(dp) :: alpha,beta
!    integer :: i, j, k, l
!
!    @:ASSERT(size(x) ==  size(A,dim = 1))
!    @:ASSERT(size(A,dim = 1) ==  size(A,dim = 2))
!    @:ASSERT(size(b) ==  size(x))
!    @:ASSERT(tol > 0.0_dp)
!
!    allocate(r(size(b),0:1))
!    allocate(rh(size(b),-1:1))
!    allocate(nu(size(b),0:1))
!    allocate(p(size(b),0:1))
!    allocate(xtmp(size(b),0:1))
!    allocate(rho(0:1))
!    allocate(omega(0:1))
!
!!    xtmp(:,0) = x
!!    r(:,0) = b(:) - matmul(A,x(:))
!!    rh(:,-1) = r(:,0)
!!    rh(:,0) = r(:,0)
!!    rho(0) = 1.0_dp
!!    alpha = 1.0_dp
!!    omega(0) = 1.0_dp
!!    nu(:,0) = 0.0_dp
!!    p(:,0) = 0.0_dp
!!    do k = 1, size(A,dim = 2)
!!      i = mod(k,2)
!!      j = mod(k+1,2)
!!      rho(i) = dot_product(rh(:,-1),rh(:,j))
!!      beta = (rho(i)/rho(j))/(alpha/omega(j))
!!      p(:,i) = r(:,j) + beta*(p(:,j)-omega(j)*nu(:,j))
!!      nu(:,i) = matmul(A,p(:,i))
!!      alpha = rho(i) / dot_product(rh(:,-1),nu(:,i))
!!      if (real(conjg(alpha)*alpha) < 1.0E-12_dp) then
!!        alpha = (1.0E-8_dp,0.0_dp)
!!      end if
!!      s(:) = r(:,j) - alpha*nu(:,i)
!!      t = matmul(A,s)
!!      omega(i) = dot_product(t,s)/dot_product(t,t)
!!      xTmp(:,i) = xTmp(:,j) + alpha * p(:,i) + omega(i) * s
!!      tmp = b(:) - matmul(A,xTmp(:,i))
!!      if (real(dot_product(tmp,tmp))<=  tol) then
!!        exit
!!      end if
!!      r(:,i) = s - omega(i) * t
!!    end do
!!
!!    x = xTmp(:,i)
!
!    do l = 1, size(A,dim = 2)
!
!      r = 0.0_dp
!      rh = 0.0_dp
!      nu = 0.0_dp
!      p = 0.0_dp
!      xtmp = 0.0_dp
!      rho = 0.0_dp
!      alpha = 0.0_dp
!      beta = 0.0_dp
!      omega = 0.0_dp
!
!      xtmp(:,0) = x
!      r(:,0) = b(:) - matmul(A,x(:))
!      rh(:,-1) = r(:,0)
!      rh(:,0) = r(:,0)
!      rho(0) = 1.0_dp
!      alpha = 1.0_dp
!      omega(0) = 1.0_dp
!      nu(:,0) = 0.0_dp
!      p(:,0) = 0.0_dp
!      do k = 1, size(A,dim = 2)
!
!        i = mod(k,2)
!        j = mod(k+1,2)
!        rho(i) = dot_product(rh(:,-1),rh(:,j))
!        beta = (rho(i)/rho(j))/(alpha/omega(j))
!        p(:,i) = r(:,j) + beta*(p(:,j)-omega(j)*nu(:,j))
!        nu(:,i) = matmul(A,p(:,i))
!        alpha = rho(i) / dot_product(rh(:,-1),nu(:,i))
!
!        s(:) = r(:,j) - alpha*nu(:,i)
!        t = matmul(A,s)
!        omega(i) = dot_product(t,s)/dot_product(t,t)
!        xTmp(:,i) = xTmp(:,j) + alpha * p(:,i) + omega(i) * s
!        tmp = b(:) - matmul(A,xTmp(:,i))
!        if (sqrt(real(dot_product(tmp,tmp))) <=  tol) then
!          exit
!        end if
!        if (real(alpha) < epsilon(1.0_dp)) exit ! about to fail on next
!        !  iteration - restart from fresh using current x values
!        r(:,i) = s - omega(i) * t
!      end do
!
!      x = xTmp(:,i)
!
!      if (sqrt(real(dot_product(tmp,tmp)))<= tol) then
!        exit
!      end if
!
!    end do
!
!  end subroutine bicgstab_cmplx


  !> Process tolerance for norm^2 of vector
  pure function setR2Tol(tol)

    !> Required tolerance, if suplied
    real(dp), intent(in), optional :: tol

    !> Tolerance
    real(dp) :: setR2Tol

    if (present(tol)) then
      setR2Tol = tol
    else
      setR2Tol = epsilon(1.0_dp)
    end if
    setR2Tol = max(setR2Tol, epsilon(1.0_dp))
    setR2Tol = setR2Tol**2

  end function setR2Tol


  !> Process an optional int value argument
  pure function setIVal(default, ii)

    !> Default
    integer, intent(in) :: default

    !> Required value, if suplied
    integer, intent(in), optional :: ii

    !> Returned value
    real(dp) :: setIVal

    setIVal = default
    if (present(ii)) then
      setiVal = ii
    end if

  end function setIVal

end module dftbp_derivs_conjgrad
