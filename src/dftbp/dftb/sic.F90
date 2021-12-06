!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines self interaction
module dftbp_dftb_sic
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_charges, only : getSummedCharges
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : mulliken, denseMulliken, getChargePerShell, getOnsitePopulation
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : addShift, totalShift
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_dftb_spin, only : getSpinShift, ud2qm, qm2ud
  use dftbp_math_blasroutines, only : symm, hemv
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: gradientMatrix ! TSicInp, TSic, TSic_init,

!
!  !> Input for the sic module
!  type TSicInp
!
!  end type TSicInp
!
!
!  !> Internal sic variables
!  type TSic
!
!  !contains
!  !  procedure ::
!
!  end type TSic

  interface gradientMatrix
    module procedure gradientMatrix_real
  end interface gradientMatrix

contains


!  !> Initializes instance.
!  subroutine TSic_init(this, inp)
!
!    !> Instance.
!    type(TSic), intent(out) :: this
!
!    !> Input data.
!    type(TSicInp), intent(in) :: inp
!
!  end subroutine TSic_init


  !> Builds anti-symmetric matrix <c_i| F^j - F^i |c_j> for spin free real case
  subroutine gradientMatrix_real(ci, env, orb, nAtom, sccCalc, overlap, sSqr, iNeighbour,&
      & nNeighbourSK, iSquare, iSparseStart, img2CentCell, species)

    !> wave function coefficients
    real(dp), intent(in) :: ci(:, :)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> number of atoms needed for atom resolved arrays
    integer, intent(in) :: nAtom

    !> Self-consistent module variables
    type(TScc), intent(inout) :: sccCalc

    !> Sparse overlap matrix
    real(dp), intent(in) :: overlap(:)

    !> Work space for dense overlap matrix
    real(dp), intent(inout) :: sSqr(:,:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    real(dp), allocatable :: M(:,:), qOrb(:,:,:), work1D(:), work2D(:,:), DeltaFock(:,:)
    real(dp), allocatable :: potAtom(:,:), potShell(:,:,:), potBlock(:,:,:,:,:), deltaPot(:,:,:,:)
    integer :: nLev, nBas, nSpin, nOrb, iLev, jLev, iSpin, iAt, iOrb1, iOrb2

    nBas = size(ci, dim=1)
    nLev = size(ci, dim=2)
    nSpin = 1

    allocate(potAtom(nAtom, nSpin))
    allocate(potShell(orb%mShell, nAtom, nSpin))
    allocate(potBlock(orb%mOrb, orb%mOrb, nAtom, nSpin, nLev))
    potBlock(:,:,:,:,:) = 0.0_dp

    call unpackHS(sSqr, overlap, iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell)

    allocate(work2D(nBas, nLev))
    allocate(qOrb(orb%mOrb,nAtom,nSpin))

    do iSpin = 1, nSpin

      call symm(work2D(:, :), 'L', sSqr, ci(:, :))

      do iLev = 1, nLev

        qOrb(:,:,:) = 0.0_dp
        do iAt = 1, nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)-1
          nOrb = iOrb2-iOrb1+1
          ! hermitian transpose used as real part only is needed
          qOrb(:nOrb,iAt,iSpin) = sum(&
              & ci(iOrb1:iOrb2, iLev) * work2D(iOrb1:iOrb2, iLev) )
        end do

        call sccCalc%updateCharges(env, qOrb, orb, species)
        call sccCalc%updateShifts(env, orb, species, iNeighbour, img2CentCell)
        potAtom(:,:) = 0.0_dp
        call sccCalc%getShiftPerAtom(potAtom(:,1))
        potShell(:,:,:) = 0.0_dp
        call sccCalc%getShiftPerL(potShell(:,:,1))

        !if (allocated(spinW)) then
        !  call getChargePerShell(qIn, orb, species, qPerShell)
        !  shellPot(:,:,:) = 0.0_dp
        !  call getSpinShift(shellPot(:,:,2:), qPerShell(:,:,2:), species, orb, spinW)
        !  potShell(:,:,2:) = potShell(:,:,2:) + shellPot(:,:,2:)
        !end if

        ! third order still to do

        call totalShift(potShell, potAtom, orb, species)
        call totalShift(potBlock(:,:,:,:,iLev), potShell, orb, species)

      end do

    end do

    deallocate(work2D)
    deallocate(qOrb)

    allocate(M(nLev, nLev))
    M(:, :) = 0.0_dp

    allocate(deltaPot(orb%mOrb, orb%mOrb, nAtom, nSpin))
    allocate(deltaFock(size(overlap),nSpin))

    allocate(work1D(nBas))

    do iLev = 1, nLev

      do jLev = 1, nLev

        deltaPot(:,:,:,:) = potBlock(:,:,:,:,iLev) - potBlock(:,:,:,:,jLev)

        deltaFock(:,:) = 0.0_dp
        call addShift(deltaFock, overlap, nNeighbourSK, iNeighbour, species, orb, iSparseStart,&
            & nAtom, img2CentCell, deltaPot)

        call unpackHS(sSqr, deltaFock(:,1), iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
            & img2CentCell)

        call hemv(work1D, sSqr, ci(:,iLev))

        M(jLev, iLev) = dot_product(ci(:, jLev), work1D)

      end do

      write(*,*)M(:, iLev)

    end do

  end subroutine gradientMatrix_real

end module dftbp_dftb_sic
