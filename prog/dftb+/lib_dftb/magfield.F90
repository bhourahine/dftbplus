!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for external magnetic fields
module HFields
  use accuracy, only : dp
  use constants
  use angmomentum, only : getLOperators
  use CommonTypes, only : TOrbitals
  implicit none

  private
  public :: shiftH, phaseHfield, unphaseHfield

  complex(dp), parameter :: im = (0.0_dp,1.0_dp)

contains


  !> Landau gauge for H along z, vector potential changing in the x direction
  function gauge(x,H)

    !> Vector potential in the gauge
    real(dp) :: gauge(3)

    !> spatial position
    real(dp), intent(in) :: x(3)

    !> magnetic field strength
    real(dp), intent(in) :: H

    gauge = 0.0_dp
    gauge(1) = -H*x(2)

  end function gauge


  !> Constructs shift potential for scalar potential part of magnetic field along z in atomic
  !> 'MKS'/'SI' units
  subroutine shiftH(shift, iShift, HFieldStrength, orb, species)

    !> block shift from the potential
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> imaginary block shift from the potential
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> magnetic field strength
    real(dp), intent(in) :: HFieldStrength

    !> momentum information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> list of the species for each atom
    integer, intent(in) :: species(:)

    integer :: iAt, nAtom, nSpin, iSp, iSh, iOrb, nSpecies
    integer :: ii, jj, kk, ll, mm, iStart, iEnd
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: speciesL(:,:,:,:)

    nAtom = size(shift,dim=3)
    nSpin = size(shift,dim=4)
    nSpecies = maxval(species(1:nAtom))

    ! spin Zeeman part
    if (nSpin > 1) then
      do iAt = 1, nAtom
        iSp = species(iAt)
        do iSh = 1, orb%nShell(iSp)
          do iOrb = orb%posShell(iSh,iSp), orb%posShell(iSh+1,iSp) -1
            shift(iOrb,iOrb,iAt,nSpin) = shift(iOrb,iOrb,iAt,nSpin)&
                & -gfac * mu_B *0.5_dp* HFieldStrength
          end do
        end do
      end do
    end if

    ! Orbital Zeeman part
    allocate(speciesL(orb%mOrb, orb%mOrb, 3, nSpecies))
    speciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))
    do ii = 1, nSpecies
      do jj = 1, orb%nShell(ii)
        Lz = 0.0_dp
        Lplus = 0.0_dp
        kk = orb%angShell(jj,ii)
        call getLOperators(kk, Lplus(:2*kk+1,:2*kk+1), Lz(:2*kk+1,:2*kk+1))
        ! x and y part not needed - H is along z
        !speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,&
        !    & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,1,ii)&
        !    & = aimag(Lplus(1:2*kk+1,1:2*kk+1))
        !speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,&
        !    & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,2,ii)&
        !    & = -real(Lplus(1:2*kk+1,1:2*kk+1))
        speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,&
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,3,ii)&
            & = aimag(Lz(1:2*kk+1,1:2*kk+1))
      end do
    end do
    deallocate(lplus)
    deallocate(lz)
    do ii = 1, natom
      isp = species(ii)
      mm = orb%norbspecies(isp)
      do jj = 1, orb%nshell(isp)
        istart = orb%posshell(jj,isp)
        iend = orb%posshell(jj+1,isp)-1
        ll = 3 ! z component only
        ishift(istart:iend,istart:iend,ii,1) =&
            & ishift(istart:iend,istart:iend,ii,1)&
            & -mu_b * hfieldstrength * speciesl(istart:iend,istart:iend,ll,isp)
      end do
    end do
    deallocate(speciesL)

  end subroutine shiftH


  !> Applies vector potential phases to transform to London orbitals for magnetic field along z in
  !> atomic 'MKS'/'SI' units
  subroutine phaseHfield(sparseRe, sparseIm, HFieldStrength, orb, coords, iNeighbor, nNeighbor,&
      & iPair, img2CentCell)

    !> coefficients of sparse matrix
    real(dp), intent(inout) :: sparseRe(:)

    !> imaginary coefficients of sparse matrix
    real(dp), intent(inout) :: sparseIm(:)

    !> magnetic field strength
    real(dp), intent(in) :: HFieldStrength

    !> Angular momentum information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> atomic positions
    real(dp), intent(in) :: coords(:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iStart, iEnd, nOrb, iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    complex(dp) :: tmp(orb%mOrb*orb%mOrb)

    nAtom = size(iNeighbor, dim=2)

  @:ASSERT(nAtom > 0)
  @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))

    tmp = 0.0_dp
    do iAtom1 = 1, nAtom
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb = orb%nOrbAtom(iAtom1)*orb%nOrbAtom(iAtom2f)
        iStart = iPair(iNeigh,iAtom1) + 1
        iEnd = iPair(iNeigh,iAtom1) + nOrb
        ! minimal coupling, CGS
        phase = exp( im * alpha_fs * 0.5_dp * dot_product(&
            & (gauge(coords(:,iAtom2),HFieldStrength)+gauge(coords(:,iAtom1),HFieldStrength)),&
            & (coords(:,iAtom2)-coords(:,iAtom1))) )
        tmp(:nOrb) = phase * (sparseRe(iStart:iEnd) + im * sparseIm(iStart:iEnd))
        sparseRe(iStart:iEnd) = real(tmp(:nOrb))
        sparseIm(iStart:iEnd) = aimag(tmp(:nOrb))
      end do
    end do

  end subroutine phaseHfield


  !> Applies vector potential phases to transform back from London orbitals
  subroutine unphaseHfield(sparseRe, sparseIm, HFieldStrength, orb, coords, iNeighbor, nNeighbor,&
      & iPair, img2CentCell)

    !> coefficients of sparse matrix
    real(dp), intent(inout) :: sparseRe(:)

    !> imaginary coefficients of sparse matrix
    real(dp), intent(inout) :: sparseIm(:)

    !> block shift from the potential
    real(dp), intent(in) :: HFieldStrength

    !> magnetic field strength
    type(TOrbitals), intent(in) :: orb

    !> Angular momentum information about the orbitals
    real(dp), intent(in) :: coords(:,:)

    !> atomic positions
    integer, intent(in) :: iNeighbor(0:,:)

    !> Neighbor list for each atom (First index from 0!)
    integer, intent(in) :: nNeighbor(:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: iPair(0:,:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: img2CentCell(:)

    !> Atomic mapping indexes.

    complex(dp) :: phase
    integer :: nAtom
    integer :: iStart, iEnd, nOrb, iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    complex(dp) :: tmp(orb%mOrb*orb%mOrb)

    nAtom = size(iNeighbor, dim=2)

  @:ASSERT(nAtom > 0)
  @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))

    tmp = 0.0_dp
    do iAtom1 = 1, nAtom
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb = orb%nOrbAtom(iAtom1)*orb%nOrbAtom(iAtom2f)
        iStart = iPair(iNeigh,iAtom1) + 1
        iEnd = iPair(iNeigh,iAtom1) + nOrb
        ! minimal coupling, CGS
        phase = exp( -im * alpha_fs * 0.5_dp * dot_product(&
            & (gauge(coords(:,iAtom2),HFieldStrength)+gauge(coords(:,iAtom1),HFieldStrength)),&
            & (coords(:,iAtom2)-coords(:,iAtom1))) )
        tmp(:nOrb) = phase * (sparseRe(iStart:iEnd) + im * sparseIm(iStart:iEnd))
        sparseRe(iStart:iEnd) = real(tmp(:nOrb))
        sparseIm(iStart:iEnd) = aimag(tmp(:nOrb))
      end do
    end do

  end subroutine unphaseHfield


end module HFields
