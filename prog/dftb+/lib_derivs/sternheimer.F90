!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using the Sternheimer equation
module dftbp_sternheimer
  use dftbp_accuracy, only : dp
  use dftbp_conjgradaxb, only : conjgradHeS_contract
  use dftbp_sparseblas, only : dftbsymm
  use dftbp_blasroutines, only : gemm, symm, herk
  implicit none

  private
  public :: sternheimerStaticReal

contains

  !> Solves -(H - e S) psi' = (H' - e S'), psi optionally with projection onto the virtual states
  subroutine sternheimerStaticReal(SP, ham, over, dHam, dOver, ei, ci, isDense, iNeighbor,&
      & nNeighbor, img2CentCell, iPair, iSquare, mOrb, dCi)

    !> Projector onto the valence states, S * P, required if states should be projected onto virtual
    !> states
    real(dp), intent(in), allocatable :: SP(:,:)

    !> Hamiltonian matrix
    real(dp), intent(in) :: ham(:)

    !> Overlap matrix
    real(dp), intent(in) :: over(:)

    !> Derivative of hamiltonian matrix
    real(dp), intent(in) :: dHam(:)

    !> Derivative of overlap matrix, if relevant to perturbation
    real(dp), intent(in), allocatable :: dOver(:)

    !> Eigenvalues
    real(dp), intent(in) :: ei(:)

    !> Eigenvectors
    real(dp), intent(in) :: ci(:,:)

    !> Should dense algebra be used
    logical, intent(in) :: isDense

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

    !> Derivative of eigenvectors
    real(dp), intent(out) :: dCi(:,:)

    integer :: ii, m, n
    real(dp), allocatable :: work(:,:), work2(:,:), ciTmp(:,:)

    m = size(ci, dim=1)
    n = size(ci, dim=2)

    allocate(ciTmp(m, n))

    if (isDense) then
      allocate(work(m,m))
    end if

    ! -H' psi
    if (isDense) then
      work(:,:) = 0.0_dp
      call unpackHS(work, dHam, iNeighbor,  nNeighbor, iSquare, iPair, img2CentCell)
      call symm(ciTmp, 'l', work, ci, alpha=-1.0_dp)
    else
      call dftbSYMM(ciTmp, dHam, ci, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,&
          & mOrb, alpha=-1.0_dp)
    end if
    if (allocated(dOver)) then
      ! add -(-e S') psi
      allocate(work2(m,n))
      do ii = 1, n
        work2(:,ii) = ei(ii) * ci(:,ii)
      end do
      if (isDense) then
        work(:,:) = 0.0_dp
        call unpackHS(work, dOver, iNeighbor,  nNeighbor, iSquare, iPair, img2CentCell)
        call symm(ciTmp, 'l', work, work2, beta=1.0_dp)
      else
        call dftbSYMM(ciTmp, dOver, work, iNeighbor, nNeighbor, img2CentCell, iPair, iSquare,&
            & mOrb, beta=1.0_dp)
      end if
      deallocate(work2)
    end if

    if (allocated(SP)) then
      allocate(work2(m,n))
      work2(:,:) = ciTmp
      call gemm(work2, SP, ciTmp, alpha = -1.0_dp, beta = 1.0_dp)
      ciTmp(:,:) = work2
      deallocate(work2)
    end if

    if (isDense) then
      work(:,:) = 0.0_dp
      call unpackHS(work, ham, iNeighbor,  nNeighbor, iSquare, iPair, img2CentCell)
      allocate(work2(m,m))
      work2(:,:) = 0.0_dp
      call unpackHS(work2, over, iNeighbor,  nNeighbor, iSquare, iPair, img2CentCell)
      call conjgradHeS_contract(work, ei, work2, ciTmp, dCi)
      deallocate(work2)
      deallocate(work)
    else
      call conjgradHeS_contract(ham, ei, over, ciTmp, dCi, iNeighbor, nNeighbor, img2CentCell,&
          & iPair, iSquare, mOrb)
    end if
    deallocate(ciTmp)

    if (allocated(SP)) then
      allocate(work(m,n))
      work(:,:) = dCi
      call gemm(work, SP, dCi, alpha = -1.0_dp, beta = 1.0_dp, transA='T')
      dCi(:,:) = work
      deallocate(work)
    end if

  end subroutine sternheimerStaticReal

end module dftbp_sternheimer
