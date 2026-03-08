!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Solve A x = b and related problems
module dftbp_math_conjugategradient
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_math_blasroutines, only : hemv, symm
  use dftbp_io_message, only : warning
  use dftbp_math_sparseblas, only : sparse_symv => symv, sparse_symm => symm
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: bicgstab, conjgrad, conjgrad_contract

  !> Conjugate gradient solution of A x = b and (H - e S)x = b
  interface conjgrad
    module procedure dense
    module procedure sparse
    module procedure dense_multix
    module procedure dense_multixHeS
    module procedure sparse_multixHeS
  end interface conjgrad

  !> Conjugate gradient solution of A x = b and (H - e S)x = b
  interface conjgrad_contract
    module procedure dense_multix_contract
    module procedure dense_multixHeS_contract
    module procedure sparse_multixHeS_contract
  end interface conjgrad_contract

  !> Biconjugate
  interface bicgstab
    module procedure bicgstab_cmplx
  end interface bicgstab

contains

  !> Calculates b in Ab=x using halting criteria |A b - x| < tol for a dense A matrix
  !! note: based on wikipedia pseudo code for conjugate gradient
  subroutine dense(A, b, x, tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> Solution vector
    real(dp), intent(in) :: b(:)

    !> Entry intial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> Halting tolerance if present (internally defaults to epsilon
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), p(size(b)), Ap(size(b))
    real(dp) :: rsold, rsnew, alpha, tol_
    integer :: ii

    @:ASSERT(size(A,dim=1) == size(A,dim=2))
    @:ASSERT(size(b) == size(x))
    @:ASSERT(size(A,dim=1) == size(x))

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    ! replaces r = b - matmul(A, x)
    r = b
    call hemv(r, A, x, 'L', alpha=-1.0_dp, beta=1.0_dp)
    p = r
    rsold = dot_product(r, r)
    do ii = 1, size(A, dim=1)
      call hemv(Ap, A, p, 'L')
      alpha = rsold / dot_product(p, Ap)
      x = x + alpha * p
      r = r - alpha * Ap
      rsnew = dot_product(r, r)
      if (rsnew < tol_) exit
      p = r + (rsnew / rsold) * p
      rsold = rsnew
    end do

  end subroutine dense


  !> Version of the above dense subroutine using the DFTB+ sparse matrix structure to store A
  subroutine sparse(A,b,x,iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb, tol)

    !> Symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> Vector
    real(dp), intent(in) :: b(:)

    !> Intial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> Neighbor list for each atom (First index from 0)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Orbitals on atoms
    type(TOrbitals), intent(in) :: orb

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), p(size(b)), Ap(size(b)), rsold, rsnew, alpha, tol_
    integer :: ii

    @:ASSERT(size(b) == size(x))

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

    ! replaces r = b - matmul(A, x)
    r = b
    call sparse_symv(r, A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb,&
        & alpha=-1.0_dp, beta=1.0_dp)
    p = r
    rsold = dot_product(r, r)
    do ii = 1, size(A, dim=1)
      ! replaces Ap = matmul(A, p)
      call sparse_symv(Ap, A, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb)
      alpha = rsold / dot_product(p, Ap)
      x = x + alpha * p
      r = r - alpha * Ap
      rsnew = dot_product(r, r)
      if (rsnew < tol_) exit
      p = r + (rsnew / rsold) * p
      rsold = rsnew
    end do

  end subroutine sparse


  !> Version of the above dense subroutine for multiple right hand sides
  !! Note: iteration of non-converged solutions are made, converged cases are simply kept unchanged
  !! (a little wasteful of GEMV calls)
  subroutine dense_multix(A, b, x, tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> Matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> On entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A,dim=1) == size(A,dim=2))
    @:ASSERT(size(A,dim=1) == n)
    @:ASSERT(all(shape(b) == [n,s]))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))

    tUnconverged = .true.
    r = b
    call symm(r, 'L', A, x, alpha=-1.0_dp, beta=1.0_dp)
    p = r
    rsold=sum(r**2,dim=1)
    do ii=1,size(A,dim=1)
      call symm(Ap, 'L', A, p)
      alpha=rsold/sum(p*Ap,dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew=sum(r**2,dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) tUnconverged(jj) = .false.
      end do
      if (all(.not.tUnconverged)) exit
      do jj = 1, s
        p(:,jj)=r(:,jj)+(rsnew(jj)/rsold(jj))*p(:,jj)
      end do
      rsold=rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)

  end subroutine dense_multix


  !> Version of the dense_multix routine above which contracts the size of the matrix of trial
  !! solution vectors as they converge during the interations
  subroutine dense_multix_contract(A,b,x,tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> Matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> On entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u
    real(dp), allocatable :: xTmp(:,:)
    integer, allocatable :: indx(:)

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(A,dim=1)==n)
    @:ASSERT(all(shape(b)==(/n,s/)))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))
    allocate(xTmp(n,s))
    allocate(indx(s))

    t = s
    xTmp = x
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged = .true.
    r = b

    call symm(r, 'L', A, xTmp, alpha=-1.0_dp, beta=1.0_dp)
    p = r
    rsold=sum(r**2,dim=1)

    do ii=1, n
      call symm(Ap(:,:t), 'L', A, p(:,:t))
      alpha(:t)=rsold(:t)/sum(p(:,:t)*Ap(:,:t),dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          xTmp(:,jj) = xTmp(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:t)=sum(r(:,:t)**2,dim=1)
      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) tUnconverged(jj) = .false.
      end do

      do jj = 1, t
        x(:,indx(jj)) = xTmp(:,jj)
      end do

      if (t > 1) then
        u = 1
        do jj = 1, t
          if (tUnconverged(jj)) then
            if (u < jj) then
              indx(u) = indx(jj)
              tUnconverged(u) = tUnconverged(jj)
              Ap(:,u) = Ap(:,jj)
              p(:,u) = p(:,jj)
              r(:,u) = r(:,jj)
              rsnew(u) = rsnew(jj)
              rsold(u) = rsold(jj)
              alpha(u) = alpha(jj)
              xTmp(:,u) = xTmp(:,jj)
            end if
            u = u +1
          end if
        end do

        t = u - 1
      end if
      if (all(.not.tUnconverged(:t)) .or. t < 1) exit

      do jj = 1, t
        p(:,jj)=r(:,jj)+(rsnew(jj)/rsold(jj))*p(:,jj)
      end do
      rsold(:t)=rsnew(:t)
    end do

    deallocate(xTmp)
    deallocate(indx)
    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)

  end subroutine dense_multix_contract


  !> Multiple right hand sides algorithm for the case (H-ES)B=X where E(i) is a parameter for each
  !! X(:,i) right hand side.
  subroutine dense_multixHeS(A,ei,over,b,x,tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> Eigenvalues
    real(dp), intent(in) :: ei(:)

    !> Overlap matrix
    real(dp), intent(in) :: over(:,:)

    !> Matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> On entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xtmp(:,:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(A,dim=1)==n)
    @:ASSERT(all(shape(A) == shape(over)))
    @:ASSERT(size(ei)==s)
    @:ASSERT(all(shape(b)==(/n,s/)))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))
    allocate(xTmp(n,s))

    tUnconverged = .true.
    r = b
    ! form Ax and over x separately
    call symm(r, 'L', A, x, alpha=-1.0_dp, beta=1.0_dp)
    call symm(xTmp, 'L', over, x)
    do jj = 1, s
      r(:,jj) = r(:,jj) + ei(jj)*xTmp(:,jj)
    end do

    p = r
    rsold=sum(r**2,dim=1)

    do ii=1,size(A,dim=1)
      call symm(Ap, 'L', A, p)
      call symm(xTmp, 'L', over, p)
      do jj = 1, s
        Ap(:,jj) = Ap(:,jj) - ei(jj)*xTmp(:,jj)
      end do
      alpha=rsold/sum(p*Ap,dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew=sum(r**2,dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) tUnconverged(jj) = .false.
      end do
      if (all(.not.tUnconverged)) exit
      do jj = 1, s
        p(:,jj)=r(:,jj)+(rsnew(jj)/rsold(jj))*p(:,jj)
      end do
      rsold=rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)
    deallocate(xTmp)

  end subroutine dense_multixHeS


  !> Version of dense_multixHeS which contracts the size of the matrix of trial solution vectors as
  !! they converge during the interations
  subroutine dense_multixHeS_contract(A,ei,over,b,x,tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> Eigenvalues
    real(dp), intent(in) :: ei(:)

    !> Overlap matrix
    real(dp), intent(in) :: over(:,:)

    !> Matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> On entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:)
    real(dp), allocatable :: rsold(:), rsnew(:), alpha(:), xStmp(:,:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u
    real(dp), allocatable :: xTmp(:,:), eiTmp(:)
    integer, allocatable :: indx(:)

    if (present(tol)) then
       tol_ = tol
    else
       tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(A,dim=1)==n)
    @:ASSERT(all(shape(A) == shape(over)))
    @:ASSERT(size(ei)==s)
    @:ASSERT(all(shape(b)==(/n,s/)))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))
    allocate(xSTmp(n,s))

    allocate(xTmp(n,s))
    allocate(indx(s))
    allocate(eiTmp(s))

    t = s
    xTmp = x
    eiTmp = ei
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged = .true.
    r = b

    call symm(r, 'L', A, x, alpha=-1.0_dp, beta=1.0_dp)
    call symm(xSTmp, 'L', over, x)
    do jj = 1, s
      r(:,jj) = r(:,jj) + eiTmp(jj)*xSTmp(:,jj)
    end do

    p = r
    rsold=sum(r**2,dim=1)

    do ii=1,4*size(A,dim=1)

      call symm(Ap(:,:t), 'L', A, p(:,:t))
      call symm(xSTmp(:,:t), 'L', over, p(:,:t))
      do jj = 1, t
        Ap(:,jj) = Ap(:,jj) - eiTmp(jj)*xSTmp(:,jj)
      end do
      alpha(:t)=rsold(:t)/sum(p(:,:t)*Ap(:,:t),dim=1)
      ! alpha=rsold/sum(p*Ap,dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          xTmp(:,jj) = xTmp(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      !rsnew=sum(r**2,dim=1)
      rsnew(:t)=sum(r(:,:t)**2,dim=1)

      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) tUnconverged(jj) = .false.
      end do

      do jj = 1, t
        x(:,indx(jj)) = xTmp(:,jj)
      end do

      if (t>1 .and. mod(ii,40)==0) then
        u = 1
        do jj = 1, t
          if (tUnconverged(jj)) then
            if (u < jj) then
              indx(u) = indx(jj)
              tUnconverged(u) = tUnconverged(jj)
              Ap(:,u) = Ap(:,jj)
              p(:,u) = p(:,jj)
              r(:,u) = r(:,jj)
              rsnew(u) = rsnew(jj)
              rsold(u) = rsold(jj)
              alpha(u) = alpha(jj)
              xTmp(:,u) = xTmp(:,jj)
              xSTmp(:,u) = xSTmp(:,jj)
              eiTmp(u) = eiTmp(jj)
            end if
            u = u +1
          end if
        end do
        !if (t /= u - 1) then
        !  write(*,*)'deflate', t, (u-1)
        !end if
        t = u - 1
      end if

      if (all(.not.tUnconverged(:t)) .or. t < 1) exit

      do jj = 1, t
        p(:,jj)=r(:,jj)+(rsnew(jj)/rsold(jj))*p(:,jj)
      end do
      rsold(:t)=rsnew(:t)

    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)
    deallocate(xSTmp)
    deallocate(xTmp)
    deallocate(eiTmp)

  end subroutine dense_multixHeS_contract


  !> Sparse version of dense_multixHeS_contract using DFTB+ matrix structure matrix of trial
  !! solution vectors as they converge during the interations
  subroutine sparse_multixHeS(A,ei,over,b,x,iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,&
      & orb, tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> Eigenvalues
    real(dp), intent(in) :: ei(:)

    !> Overlap matrix
    real(dp), intent(in) :: over(:)

    !> Matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> On entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Maximal number of orbitals on an atom.
    type(TOrbitals), intent(in) :: orb

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xtmp(:,:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(all(shape(A) == shape(over)))
    @:ASSERT(size(ei)==s)
    @:ASSERT(all(shape(b)==(/n,s/)))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))
    allocate(xTmp(n,s))

    tUnconverged = .true.
    r = b
    call sparse_symm(r, 'L', A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb,&
        & alpha=-1.0_dp, beta=1.0_dp)
    ! call symm(r,A,x,alpha=-1.0_dp,beta=1.0_dp)
    call sparse_symm(xTmp, 'L', over, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb)
    ! call symm(xTmp,over,x)
    do jj = 1, s
      r(:,jj) = r(:,jj) + ei(jj)*xTmp(:,jj)
    end do

    p = r
    rsold=sum(r**2,dim=1)

    do ii=1,4*size(ei)
      call sparse_symm(Ap, 'L', A, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb)
      call sparse_symm(xTmp, 'L', over, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb)
      do jj = 1, s
        Ap(:,jj) = Ap(:,jj) - ei(jj)*xTmp(:,jj)
      end do
      alpha=rsold/sum(p*Ap,dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew=sum(r**2,dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) tUnconverged(jj) = .false.
      end do
      if (all(.not.tUnconverged)) exit
      do jj = 1, s
        p(:,jj)=r(:,jj)+(rsnew(jj)/rsold(jj))*p(:,jj)
      end do
      rsold=rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)
    deallocate(xTmp)

  end subroutine sparse_multixHeS


  !> Sparse version of dense_multixHeS_contract using DFTB+ matrix structure matrix of trial
  !! solution vectors as they converge during the interations
  subroutine sparse_multixHeS_contract(A,ei,over,b,x,iNeighbor, nNeighbor, img2CentCell, iPair,&
      & iSquare, orb, tol)

    !> Square symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> Eigenvalues
    real(dp), intent(in) :: ei(:)

    !> Overlap matrix
    real(dp), intent(in) :: over(:)

    !> Matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> On entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol
    
    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:)
    real(dp), allocatable :: rsold(:), rsnew(:), alpha(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u
    real(dp), allocatable :: xTmp(:,:)
    integer, allocatable :: indx(:)
    character(lc) :: lcTmp

    if (present(tol)) then
       tol_ = tol
    else
       tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(all(shape(b)==(/n,s/)))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))

    allocate(xTmp(n,s))
    allocate(indx(s))

    t = s
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged = .true.

    r = b
    call sparse_symm(r, 'L', A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb,&
        & alpha=-1.0_dp, beta=1.0_dp)
    call sparse_symm(xtmp, 'L', over,x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, orb)
    ! call symm(r,A,x,alpha=-1.0_dp,beta=1.0_dp)
    ! call symm(xtmp,over,x)
    do jj = 1, t
      r(:,jj) = r(:,jj) + ei(indx(jj))*xTmp(:,jj)
    end do

    p = r
    rsold=sum(r**2,dim=1)

    do ii = 1, 4 * size(ei)
      !if (mod(ii,5)==0) then
      !  write(*,*)'CG cycle',ii,'of',4*size(ei)
      !  write(*,*)'Converged',count(.not.tUnconverged(:t)),'/',t
      !end if
      call sparse_symm(Ap(:,:t), 'L', A, p(:,:t), iNeighbor, nNeighbor, img2CentCell, iPair,&
          & iSquare, orb)
      call sparse_symm(xtmp(:,:t), 'L', over, p(:,:t), iNeighbor, nNeighbor, img2CentCell, iPair,&
          & iSquare, orb)
      ! call symm(Ap(:,:t),A,p(:,:t))
      ! call symm(xtmp(:,:t),over,p(:,:t))
      do jj = 1, t
        Ap(:,jj) = Ap(:,jj) - ei(indx(jj))*xTmp(:,jj)
      end do
      alpha(:t)=rsold(:t)/sum(p(:,:t)*Ap(:,:t),dim=1)
      ! alpha=rsold/sum(p*Ap,dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          x(:,indx(jj)) = x(:,indx(jj)) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      !rsnew=sum(r**2,dim=1)
      rsnew(:t)=sum(r(:,:t)**2,dim=1)

      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) tUnconverged(jj) = .false.
      end do

      if (t>1 .and. mod(ii,40)==0) then
        u = 1
        do jj = 1, t
          if (tUnconverged(jj)) then
            if (u < jj) then
              indx(u) = indx(jj)
              tUnconverged(u) = tUnconverged(jj)
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
        !if (t /= u - 1) then
        !  write(*,*)'deflate', t, '=>', (u-1)
        !end if
        t = u - 1
      end if

      if (all(.not.tUnconverged(:t)) .or. t < 1) exit

      do jj = 1, t
        p(:,jj)=r(:,jj)+(rsnew(jj)/rsold(jj))*p(:,jj)
      end do
      rsold(:t)=rsnew(:t)

    end do


    if (any(tUnconverged)) then
      write(lcTmp,"('Not all vectors converged in CG ',I0, ' failed')")count(tUnconverged)
      call warning(trim(lcTmp))
    end if

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)
    deallocate(xTmp)

  end subroutine sparse_multixHeS_contract


  !> Modified version of BICGSTAB algorithm to solve Ax=b, where the interation restarts itself if
  !! alpha goes to 0 (so if it were about to fail on the next iteration).
  subroutine bicgstab_cmplx(a,b,x,tol)

    !> On entry intial guess to solution, on exit the solution
    complex(dp), intent(in) :: A(:,:)

    !> Matrix
    complex(dp), intent(in) :: b(:)

    !> Right hand side
    complex(dp), intent(inout) :: x(:)

    !> Tolerance for | Ax - b |
    real(dp), intent(in) :: tol

    complex(dp), allocatable :: r(:,:), rh(:,:), nu(:,:), p(:,:), xtmp(:,:), rho(:), omega(:)
    complex(dp) :: s(size(b)), t(size(b)), tmp(size(b)), alpha,beta
    integer :: ii, jj, kk, ll

    @:ASSERT(size(x)==size(A,dim=1))
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(b)==size(x))
    @:ASSERT(tol > 0.0_dp)

    allocate(r(size(b),0:1))
    allocate(rh(size(b),-1:1))
    allocate(nu(size(b),0:1))
    allocate(p(size(b),0:1))
    allocate(xtmp(size(b),0:1))
    allocate(rho(0:1))
    allocate(omega(0:1))

    do ll = 1, size(A,dim=2)

      r = 0.0_dp
      rh = 0.0_dp
      nu = 0.0_dp
      p = 0.0_dp
      xtmp = 0.0_dp
      rho = 0.0_dp
      alpha = 0.0_dp
      beta = 0.0_dp
      omega = 0.0_dp

      xtmp(:,0) = x
      r(:,0) = b(:) - matmul(A,x(:))
      rh(:,-1) = r(:,0)
      rh(:,0) = r(:,0)
      rho(0) = 1.0_dp
      alpha = 1.0_dp
      omega(0) = 1.0_dp
      nu(:,0) = 0.0_dp
      p(:,0) = 0.0_dp
      do kk = 1, size(A,dim=2)

        ii = mod(kk,2)
        jj = mod(kk+1,2)
        rho(ii) = dot_product(rh(:,-1),rh(:,jj))
        beta = (rho(ii)/rho(jj))/(alpha/omega(jj))
        p(:,ii) = r(:,jj) + beta*(p(:,jj)-omega(jj)*nu(:,jj))
        nu(:,ii) = matmul(A,p(:,ii))
        alpha = rho(ii) / dot_product(rh(:,-1),nu(:,ii))

        s(:) = r(:,jj) - alpha*nu(:,ii)
        t = matmul(A,s)
        omega(ii) = dot_product(t,s)/dot_product(t,t)
        xTmp(:,ii) = xTmp(:,jj) + alpha * p(:,ii) + omega(ii) * s
        tmp = b(:) - matmul(A,xTmp(:,ii))
        if (sqrt(real(dot_product(tmp,tmp))) <= tol) exit
        if (real(alpha) < epsilon(1.0_dp)) exit ! about to fail on next
        !  iteration - restart from fresh using current x values
        r(:,ii) = s - omega(ii) * t
      end do

      x = xTmp(:,ii)

      if (sqrt(real(dot_product(tmp,tmp)))<= tol) exit

    end do

    deallocate(r)
    deallocate(rh)
    deallocate(nu)
    deallocate(p)
    deallocate(xtmp)
    deallocate(rho)
    deallocate(omega)

  end subroutine bicgstab_cmplx
  
end module dftbp_math_conjugategradient
