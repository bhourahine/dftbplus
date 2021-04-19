!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Conjugate gradient solvers for problems of the forms A x = b
module dftbp_conjgradaxb
  use dftbp_accuracy, only : dp, lc
  use dftbp_assert
  use dftbp_blasroutines, only : gemm, symm
  use dftbp_message, only : warning
  use dftbp_sparseblas, only : dftbsymm, dftbsymv
  implicit none

  private
  public :: conjgrad, conjgrad_contract, conjgradHeS, conjgradHeS_contract, bicgstab

  !> Conjugate gradient
  interface conjgrad
    module procedure dense
    module procedure sparse
    module procedure dense_multiRHS
    module procedure sparse_multiRHS
    module procedure densePrecond
  end interface conjgrad

  !> Conjugate gradient for multiple RHS, contracting problem as different vectors converge
  interface conjgrad_contract
    module procedure dense_multiRHS_contract
    module procedure sparse_multiRHS_contract
  end interface conjgrad_contract

  !> Solves for (H - e_i S) x_i = b_i for multiple RHS
  interface conjgradHeS
    module procedure dense_multiRHS_HeS
    module procedure sparse_multiRHS_HeS
  end interface conjgradHeS

  !> Solves for (H - e_i S) x_i = b_i for multiple RHS, contracting problem as different vectors
  !> converge
  interface conjgradHeS_contract
    module procedure dense_multiRHS_HeS_contract
    module procedure sparse_multiRHS_HeS_contract
  end interface conjgradHeS_contract

  !> Biconjugate (stabilized) gradient method
  interface bicgstab
    module procedure bicgstab_cmplx
  end interface bicgstab

contains

  !> Calculates x in A x = b using halting criteria |A x - b| < tol
  !> note based on wikipedia pseudocode
  subroutine dense(A, b, x, tol)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> RHS
    real(dp), intent(in) :: b(:)

    !> on entry initial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), p(size(b)), Ap(size(b)), rsold, rsnew, alpha, tol_
    integer :: ii

    @:ASSERT(size(A,dim=1) == size(A,dim=2))
    @:ASSERT(size(b) == size(x))
    @:ASSERT(size(A,dim=1) == size(x))

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

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
      p(:) = r + (rsnew / rsold) * p ! Fletcher-Reeves
      rsold = rsnew
    end do

  end subroutine dense


  !> Version of the dense subroutine using the DFTB+ sparse matrix structure to store A
  subroutine sparse(A, b, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb, tol)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> RHS
    real(dp), intent(in) :: b(:)

    !> on entry intial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), p(size(b)), Ap(size(b)), rsold, rsnew, alpha, tol_
    integer :: ii

    @:ASSERT(size(b)==size(x))

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

    ! replaces r=b-matmul(A,x)
    r(:) = b
    call dftbSYMV(r, A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb, alpha=-1.0_dp,&
        & beta=1.0_dp)
    p(:) = r
    rsold = dot_product(r, r)
    do ii = 1, size(A, dim=1)
      ! replaces Ap=matmul(A,p)
      call dftbSYMV(Ap, A, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb)
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

  end subroutine sparse


  !> Version of the dense subroutine for multiple right hand sides
  !>
  !> note Iteration of non-converged solutions are made, converged cases are simply kept unchanged
  !> (a little wasteful of SYMM calls)
  subroutine dense_multiRHS(A, b, x, tol)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> matrix of RHS vectors
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
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
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A, dim=1) == size(A, dim=2))
    @:ASSERT(size(A, dim=1) == n)
    @:ASSERT(all(shape(b) == [n,s]))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))

    tUnconverged(:) = .true.
    r(:,:) = b
    call symm(r, 'l', A, x, alpha=-1.0_dp, beta=1.0_dp)
    p(:,:) = r
    rsold = sum(r**2, dim=1)
    do ii = 1, size(A, dim=1)
      call symm(Ap, 'l', A, p)
      alpha(:) = rsold / sum(p * Ap, dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:) = sum(r**2, dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj) < tol_) then
          tUnconverged(jj) = .false.
        end if
      end do
      if (all(.not.tUnconverged)) then
        exit
      end if
      do jj = 1, s
        p(:,jj) = r(:,jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:) = rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)

  end subroutine dense_multiRHS


  !> Version of the sparse subroutine for multiple right hand sides
  !>
  !> note Iteration of non-converged solutions are made, converged cases are simply kept unchanged
  !> (a little wasteful of GEMV calls)
  subroutine sparse_multiRHS(A, b, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb, tol)

    !> square symmetric positive definite matrix in sparse storage
    real(dp), intent(in) :: A(:)

    !> matrix of RHS vectors
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> halting tolerance if present (internally defaults to epsilon otherwise)
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
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(all(shape(b) == [n, s]))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n,s))
    allocate(tUnconverged(s))

    tUnconverged(:) = .true.
    r(:,:) = b
    call dftbSYMM(r, A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb, alpha=-1.0_dp,&
        & beta=1.0_dp)
    p(:,:) = r
    rsold = sum(r**2, dim=1)
    do ii = 1, size(A, dim=1)
      call dftbSYMM(Ap, A, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
      alpha(:) = rsold / sum(p * Ap, dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:) = sum(r**2, dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj) < tol_) then
          tUnconverged(jj) = .false.
        end if
      end do
      if (all(.not.tUnconverged)) then
        exit
      end if
      do jj = 1, s
        p(:,jj) = r(:,jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:) = rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)

  end subroutine sparse_multiRHS


  !> Multiple right hand sides algorithm for the case (H-ES) X = B where E(i) is a parameter for
  !> each X(:,i) LHS.
  subroutine dense_multiRHS_HeS(A,ei,over,b,x,tol)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> E vector of values
    real(dp), intent(in) :: ei(:)

    !> S matrix
    real(dp), intent(in) :: over(:,:)

    !> RHS matrix
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit right hand side solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
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
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

    n = size(x, dim=1)
    s = size(x, dim=2)
    @:ASSERT(size(A,dim=1) == size(A,dim=2))
    @:ASSERT(size(A,dim=1) == n)
    @:ASSERT(all(shape(A) == shape(over)))
    @:ASSERT(size(ei) == s)
    @:ASSERT(all(shape(b) == [n,s]))

    allocate(alpha(s))
    allocate(r(n,s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n,s))
    allocate(Ap(n, s))
    allocate(tUnconverged(s))
    allocate(xTmp(n, s))

    tUnconverged(:) = .true.
    r(:,:) = b
    ! form Ax and over x separately
    call symm(r, 'l', A, x, alpha=-1.0_dp, beta=1.0_dp)
    call symm(xTmp, 'l', over, x)
    do jj = 1, s
      r(:,jj) = r(:,jj) + ei(jj) * xTmp(:,jj)
    end do

    p(:,:) = r
    rsold=sum(r**2, dim=1)

    do ii = 1, size(A, dim=1)
      call symm(Ap, 'l', A, p)
      call symm(xTmp, 'l', over, p)
      do jj = 1, s
        Ap(:,jj) = Ap(:,jj) - ei(jj) * xTmp(:,jj)
      end do
      alpha = rsold / sum(p * Ap, dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew = sum(r**2, dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) then
          tUnconverged(jj) = .false.
        end if
      end do
      if (all(.not.tUnconverged)) exit
      do jj = 1, s
        p(:,jj) = r(:,jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:) = rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)
    deallocate(xTmp)

  end subroutine dense_multiRHS_HeS


  !> Sparse version of dense_multiRHS_HeS using DFTB+ matrix structure matrix of trial solution
  !> vectors as they converge during the interations
  subroutine sparse_multiRHS_HeS(A,ei,over,b,x,iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,&
      & mOrb, tol)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> E vector of values
    real(dp), intent(in) :: ei(:)

    !> S matrix
    real(dp), intent(in) :: over(:)

    !> matrix of RHS
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> halting tolerance if present (internally defaults to epsilon otherwise)
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
    @:ASSERT(size(ei) == s)
    @:ASSERT(all(shape(b) == [n,s]))

    allocate(alpha(s))
    allocate(r(n, s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n, s))
    allocate(Ap(n, s))
    allocate(tUnconverged(s))
    allocate(xTmp(n, s))

    tUnconverged(:) = .true.
    r(:,:) = b
    call dftbSYMM(r, A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb, alpha=-1.0_dp,&
        & beta=1.0_dp)
    call dftbSYMM(xTmp, over, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb)
    do jj = 1, s
      r(:,jj) = r(:,jj) + ei(jj) * xTmp(:,jj)
    end do

    p(:,:) = r
    rsold = sum(r**2, dim=1)

    do ii=1,4*size(ei)
      call dftbSYMM(Ap, A, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
      call dftbSYMM(xTmp, over, p, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
      do jj = 1, s
        Ap(:,jj) = Ap(:,jj) - ei(jj) * xTmp(:,jj)
      end do
      alpha = rsold / sum(p * Ap, dim=1)
      do jj = 1, s
        if (tUnconverged(jj)) then
          x(:,jj) = x(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:) = sum(r**2, dim=1)
      do jj = 1, s
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) then
          tUnconverged(jj) = .false.
        end if
      end do
      if (all(.not.tUnconverged)) exit
      do jj = 1, s
        p(:, jj) = r(:, jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:) = rsnew
    end do

    deallocate(Ap)
    deallocate(rsold)
    deallocate(rsnew)
    deallocate(p)
    deallocate(r)
    deallocate(alpha)
    deallocate(tUnconverged)
    deallocate(xTmp)

  end subroutine sparse_multiRHS_HeS


  !> Calculates x in A x = b with preconditioning and using halting criteria |A x - b| < tol
  !> note: based on wikipedia pseudocode
  subroutine densePrecond(A, b, x, Minv, tol)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> RHS
    real(dp), intent(in) :: b(:)

    !> on entry initial guess to solution, on exit the solution
    real(dp), intent(inout) :: x(:)

    !> Preconditioning matrix for problem (square)
    real(dp), intent(in) :: Minv(:,:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    real(dp) :: r(size(b)), z(size(b)), zOld(size(b)), p(size(b)), Ap(size(b)), alpha, beta, tol_
    integer :: ii

    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(b)==size(x))
    @:ASSERT(size(A,dim=1)==size(x))
    @:ASSERT(all(shape(A) == shape(Minv)))

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_, epsilon(1.0_dp))
    tol_ = tol_**2

    r(:) = b - matmul(A, x)
    zOld(:) = matmul(Minv, r)
    p(:) = zOld
    do ii = 1, size(A, dim=1)
      Ap(:) = matmul(A, p)
      alpha = dot_product(r,zOld) / dot_product(p, Ap)
      x(:) = x + alpha * p
      r(:) = r - alpha * Ap
      if (dot_product(r, r) < tol_) then
        exit
      end if
      z(:) = matmul(Minv, r)
      beta = dot_product(r, z - zOld) / dot_product(r, zOld) ! Polak-Ribière
      p(:) = z + beta * p
      zOld(:) = z
    end do

  end subroutine densePrecond


  !> Version of the dense_multiRHS routine above which contracts the size of the matrix of trial
  !> solution vectors as they converge during the interations
  subroutine dense_multiRHS_contract(A, b, x, tol, nIterContract)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> matrix of RHS vectors
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    !> Period at which to check if some vectors have converged (and reduce size of problem)
    integer, intent(in), optional :: nIterContract

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xTmp(:,:)
    integer, allocatable :: indx(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u, nConvIt

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    if (present(nIterContract)) then
      nConvIt = nIterContract
    else
      nConvIt = 40
    end if

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(A,dim=1)==n)
    @:ASSERT(all(shape(b)==[n,s]))

    allocate(alpha(s))
    allocate(r(n, s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n, s))
    allocate(Ap(n, s))
    allocate(tUnconverged(s))
    allocate(xTmp(n, s))
    allocate(indx(s))

    t = s
    xTmp(:,:) = x
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged(:) = .true.
    r(:,:) = b

    call symm(r, 'l', A, xTmp, alpha=-1.0_dp, beta=1.0_dp)
    p(:,:) = r
    rsold(:) = sum(r**2,dim=1)

    do ii = 1, n
      call symm(Ap(:,:t), 'l', A, p(:, :t))
      alpha(:t) = rsold(:t) / sum(p(:, :t) * Ap(:, :t), dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          xTmp(:,jj) = xTmp(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:t) = sum(r(:, :t)**2, dim=1)
      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj) < tol_) then
          tUnconverged(jj) = .false.
        end if
      end do

      do jj = 1, t
        x(:, indx(jj)) = xTmp(:, jj)
      end do

      if (t > 1 .and. mod(ii, nConvIt)==0) then
        ! remove converged LHS vectors
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
        p(:,jj) = r(:,jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:t) = rsnew(:t)
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

  end subroutine dense_multiRHS_contract


  !> Version of the sparse_multiRHS routine above which contracts the size of the matrix of trial
  !> solution vectors as they converge during the interations
  subroutine sparse_multiRHS_contract(A, b, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,&
      & mOrb, tol, nIterContract)

    !> sparse symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> matrix of RHS vectors
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    !> Period at which to check if some vectors have converged (and reduce size of problem)
    integer, intent(in), optional :: nIterContract

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xTmp(:,:)
    integer, allocatable :: indx(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u, nConvIt

    if (present(tol)) then
      tol_ = tol
    else
      tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    if (present(nIterContract)) then
      nConvIt = nIterContract
    else
      nConvIt = 40
    end if

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(all(shape(b)==[n,s]))

    allocate(alpha(s))
    allocate(r(n, s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n, s))
    allocate(Ap(n, s))
    allocate(tUnconverged(s))
    allocate(xTmp(n, s))
    allocate(indx(s))

    t = s
    xTmp(:,:) = x
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged(:) = .true.
    r(:,:) = b
    call dftbSYMM(r, A, xTmp, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb,&
        & alpha=-1.0_dp, beta=1.0_dp)
    p(:,:) = r
    rsold(:) = sum(r**2,dim=1)

    do ii = 1, n
      call dftbSYMM(Ap(:,:t), A, p(:, :t), iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
      alpha(:t) = rsold(:t) / sum(p(:, :t) * Ap(:, :t), dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          xTmp(:,jj) = xTmp(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:t) = sum(r(:, :t)**2, dim=1)
      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj) < tol_) then
          tUnconverged(jj) = .false.
        end if
      end do

      do jj = 1, t
        x(:, indx(jj)) = xTmp(:, jj)
      end do

      if (t > 1 .and. mod(ii, nConvIt)==0) then
        ! remove converged LHS vectors
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
        p(:,jj) = r(:,jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:t) = rsnew(:t)
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

  end subroutine sparse_multiRHS_contract


  !> Version of dense_multiRHS_HeS which contracts the size of the matrix of trial solution vectors
  !> as they converge during the interations
  subroutine dense_multiRHS_HeS_contract(A, ei, over, b, x, tol, nIterContract)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:,:)

    !> E vector of values
    real(dp), intent(in) :: ei(:)

    !> S matrix
    real(dp), intent(in) :: over(:,:)

    !> matrix of solution vectors
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    !> Period at which to check if some vectors have converged (and reduce size of problem)
    integer, intent(in), optional :: nIterContract

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xStmp(:,:)
    real(dp), allocatable :: xTmp(:,:), eiTmp(:)
    integer, allocatable :: indx(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u, nConvIt
    character(lc) :: lcTmp

    if (present(tol)) then
       tol_ = tol
    else
       tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2
    if (present(nIterContract)) then
      nConvIt = nIterContract
    else
      nConvIt = 40
    end if

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(A,dim=1)==n)
    @:ASSERT(all(shape(A) == shape(over)))
    @:ASSERT(size(ei)==s)
    @:ASSERT(all(shape(b)==[n,s]))

    allocate(alpha(s))
    allocate(r(n, s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n, s))
    allocate(Ap(n, s))
    allocate(tUnconverged(s))
    allocate(xSTmp(n, s))

    allocate(xTmp(n, s))
    allocate(indx(s))
    allocate(eiTmp(s))

    t = s
    xTmp(:,:) = x
    eiTmp(:) = ei
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged(:) = .true.
    r(:,:) = b

    call symm(r, 'l', A, x, alpha=-1.0_dp, beta=1.0_dp)
    call symm(xSTmp, 'l', over, x)
    do jj = 1, s
      r(:,jj) = r(:, jj) + eiTmp(jj) * xSTmp(:, jj)
    end do

    p(:,:) = r
    rsold(:) = sum(r**2, dim=1)

    do ii = 1, 4*size(A, dim=1)

      call symm(Ap(:, :t), 'l', A, p(:, :t))
      call symm(xSTmp(:, :t), 'l', over, p(:, :t))
      do jj = 1, t
        Ap(:,jj) = Ap(:,jj) - eiTmp(jj)*xSTmp(:,jj)
      end do
      alpha(:t) = rsold(:t) / sum(p(:, :t) * Ap(:, :t), dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          xTmp(:,jj) = xTmp(:,jj) + alpha(jj) * p(:,jj)
          r(:,jj) = r(:,jj) - alpha(jj) * Ap(:,jj)
        end if
      end do
      rsnew(:t) = sum(r(:, :t)**2, dim=1)

      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) then
          tUnconverged(jj) = .false.
        end if
      end do

      do jj = 1, t
        x(:, indx(jj)) = xTmp(:, jj)
      end do

      if (t > 1 .and. mod(ii, nConvIt) == 0) then
        ! remove converged LHS vectors
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
        t = u - 1
      end if

      if (all(.not.tUnconverged(:t)) .or. t < 1) then
        exit
      end if

      do jj = 1, t
        p(:, jj) = r(:, jj) + (rsnew(jj) / rsold(jj)) * p(:, jj)
      end do
      rsold(:t) = rsnew(:t)

    end do

    if (any(tUnconverged)) then
      write(lcTmp,"('Not all vectors converged in CG: ', I0, ' failed')")count(tUnconverged)
      call warning(trim(lcTmp))
    end if

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

  end subroutine dense_multiRHS_HeS_contract


  !> Sparse version of dense_multiRHS_HeS_contract using DFTB+ sparse matrix. Contracting
  !> trial solution vectors as they converge during the interations
  subroutine sparse_multiRHS_HeS_contract(A,ei,over,b,x,iNeighbor, nNeighbor, img2CentCell, iPair,&
      & iSquare, mOrb, tol, nIterContract)

    !> square symmetric positive definite matrix
    real(dp), intent(in) :: A(:)

    !> E vector of values
    real(dp), intent(in) :: ei(:)

    !> S matrix
    real(dp), intent(in) :: over(:)

    !> RHS
    real(dp), intent(in) :: b(:,:)

    !> on entry matrix of initial guesses for solutions, on exit solution vectors
    real(dp), intent(inout) :: x(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> halting tolerance if present (internally defaults to epsilon otherwise)
    real(dp), intent(in), optional :: tol

    !> Period at which to check if some vectors have converged (and reduce size of problem)
    integer, intent(in), optional :: nIterContract

    real(dp), allocatable :: r(:,:), p(:,:), Ap(:,:), rsold(:), rsnew(:), alpha(:), xTmp(:,:)
    integer, allocatable :: indx(:)
    logical, allocatable :: tUnconverged(:)
    real(dp) :: tol_
    integer :: ii, jj, s, n, t, u, nConvIt

    character(lc) :: lcTmp

    if (present(tol)) then
       tol_ = tol
    else
       tol_ = epsilon(1.0_dp)
    end if
    tol_ = max(tol_,epsilon(1.0_dp))
    tol_ = tol_**2

    if (present(nIterContract)) then
      nConvIt = nIterContract
    else
      nConvIt = 40
    end if

    n = size(x,dim=1)
    s = size(x,dim=2)
    @:ASSERT(all(shape(b)==[n,s]))

    allocate(alpha(s))
    allocate(r(n, s))
    allocate(rsold(s))
    allocate(rsnew(s))
    allocate(p(n, s))
    allocate(Ap(n, s))
    allocate(tUnconverged(s))

    allocate(xTmp(n, s))
    allocate(indx(s))

    t = s
    do ii = 1, t
      indx(ii) = ii
    end do

    tUnconverged(:) = .true.

    r(:,:) = b
    call dftbSYMM(r, A, x, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare, mOrb, alpha=-1.0_dp,&
        & beta=1.0_dp)
    call dftbSYMM(xtmp,over,x,iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
    do jj = 1, t
      r(:,jj) = r(:,jj) + ei(indx(jj))*xTmp(:,jj)
    end do

    p(:,:) = r
    rsold(:) = sum(r**2, dim=1)

    do ii = 1, 4*size(ei)
      call dftbSYMM(Ap(:,:t),A,p(:,:t),iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
      call dftbSYMM(xtmp(:,:t),over,p(:,:t),iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,mOrb)
      do jj = 1, t
        Ap(:,jj) = Ap(:,jj) - ei(indx(jj)) * xTmp(:,jj)
      end do
      alpha(:t) = rsold(:t) / sum(p(:, :t) * Ap(:, :t), dim=1)
      do jj = 1, t
        if (tUnconverged(jj)) then
          x(:, indx(jj)) = x(:, indx(jj)) + alpha(jj) * p(:, jj)
          r(:, jj) = r(:, jj) - alpha(jj) * Ap(:, jj)
        end if
      end do
      rsnew(:t) = sum(r(:, :t)**2, dim=1)

      do jj = 1, t
        if (tUnconverged(jj) .and. rsnew(jj)<tol_) then
          tUnconverged(jj) = .false.
        end if
      end do

      if (t>1 .and. mod(ii, nConvIt)==0) then
        ! remove converged LHS vectors
        u = 1
        do jj = 1, t
          if (tUnconverged(jj)) then
            if (u < jj) then
              indx(u) = indx(jj)
              tUnconverged(u) = tUnconverged(jj)
              Ap(:, u) = Ap(:, jj)
              p(:, u) = p(:, jj)
              r(:, u) = r(:, jj)
              rsnew(u) = rsnew(jj)
              rsold(u) = rsold(jj)
              alpha(u) = alpha(jj)
            end if
            u = u +1
          end if
        end do
        t = u - 1
      end if

      if (all(.not.tUnconverged(:t)) .or. t < 1) then
        exit
      end if

      do jj = 1, t
        p(:, jj) = r(:, jj) + (rsnew(jj) / rsold(jj)) * p(:,jj)
      end do
      rsold(:t) = rsnew(:t)

    end do

    if (any(tUnconverged)) then
      write(lcTmp,"('Not all vectors converged in CG: ', I0, ' failed')")count(tUnconverged)
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

  end subroutine sparse_multiRHS_HeS_contract


  !> Modified version of BICGSTAB algorithm to solve Ax=b, where the interation restarts itself if
  !> alpha goes to 0 (which happens if it were about to fail on the next iteration).
  subroutine bicgstab_cmplx(a, b, x, tol)

    !> Matrix
    complex(dp), intent(in) :: A(:,:)

    !> RHS
    complex(dp), intent(in) :: b(:)

    !> on entry intial guess to solution, on exit the solution
    complex(dp), intent(inout) :: x(:)

    !> tolerance for | Ax - b |
    real(dp), intent(in) :: tol

    complex(dp), allocatable :: r(:,:), rh(:,:), nu(:,:), p(:,:), xtmp(:,:)
    complex(dp), allocatable :: rho(:), omega(:)
    complex(dp) :: s(size(b)), t(size(b)), tmp(size(b))
    complex(dp) :: alpha, beta
    integer :: i, j, k, l
    real(dp) :: tol2

    @:ASSERT(size(x)==size(A,dim=1))
    @:ASSERT(size(A,dim=1)==size(A,dim=2))
    @:ASSERT(size(b)==size(x))
    @:ASSERT(tol > 0.0_dp)

    tol2 = tol**2

    allocate(r(size(b), 0:1))
    allocate(rh(size(b), -1:1))
    allocate(nu(size(b), 0:1))
    allocate(p(size(b), 0:1))
    allocate(xtmp(size(b), 0:1))
    allocate(rho(0:1))
    allocate(omega(0:1))

    do l = 1, size(A, dim=2)

      r(:,:) = 0.0_dp
      rh(:,:) = 0.0_dp
      nu(:,:) = 0.0_dp
      p(:,:) = 0.0_dp
      xtmp(:,:) = 0.0_dp
      rho(:) = 0.0_dp
      alpha = 0.0_dp
      beta = 0.0_dp
      omega(:) = 0.0_dp

      xtmp(:,0) = x
      r(:,0) = b(:) - matmul(A,x(:))
      rh(:,-1) = r(:,0)
      rh(:,0) = r(:,0)
      rho(0) = 1.0_dp
      alpha = 1.0_dp
      omega(0) = 1.0_dp
      nu(:,0) = 0.0_dp
      p(:,0) = 0.0_dp

      do k = 1, size(A, dim=2)
        i = mod(k, 2)
        j = mod(k+1, 2)
        rho(i) = dot_product(rh(:, -1), rh(:, j))
        beta = (rho(i) / rho(j)) / (alpha / omega(j))
        p(:,i) = r(:, j) + beta * (p(:, j) - omega(j) * nu(:, j))
        nu(:,i) = matmul(A, p(:, i))
        alpha = rho(i) / dot_product(rh(:, -1), nu(:, i))

        s(:) = r(:, j) - alpha * nu(:, i)
        t = matmul(A, s)
        omega(i) = dot_product(t, s) / dot_product(t, t)
        xTmp(:,i) = xTmp(:, j) + alpha * p(:, i) + omega(i) * s
        tmp = b(:) - matmul(A, xTmp(:, i))
        if (real(dot_product(tmp, tmp)) <= tol2) then
          exit
        end if
        if (real(alpha) < epsilon(1.0_dp)) then
          ! about to fail on next iteration -> restart from fresh using current x values
          exit
        end if
        r(:, i) = s - omega(i) * t
      end do

      x = xTmp(:, i)

      if (real(dot_product(tmp, tmp)) <= tol2) then
        exit
      end if

    end do

    deallocate(r)
    deallocate(rh)
    deallocate(nu)
    deallocate(p)
    deallocate(xtmp)
    deallocate(rho)
    deallocate(omega)

  end subroutine bicgstab_cmplx

end module dftbp_conjgradaxb
