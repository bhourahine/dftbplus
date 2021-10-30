!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for external magnetic fields - static magnetic field
module dftbp_dftb_magfield
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : gfac, mu_B, pi, hbar, Coulomb__au, AA__Bohr
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_math_simplealgebra, only : cross3
  use dftbp_dftb_vectorpots, only : circular, xyPlane
  use dftbp_math_angmomentum, only : getLOperatorsForSpecies
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_integral, only : TIntegral
  implicit none

  private
  public :: TMagField, TMagField_init, TMagFieldInp


  !> Input for magnetic field
  type TMagFieldInp

    !> If periodic, number of unit cell repeats along first lattice vector such that one flux quanta
    !> passes through the resulting supercell
    real(dp), allocatable :: nQuanta

    !> Magnetic field as vector (not periodic system)
    real(dp), allocatable :: B(:)

    !> Origin for circular gauge magnetic field (not periodic system)
    real(dp), allocatable :: origin(:)

    !> Are only the scalar (gauge independent) effects included
    logical :: isVectorPotentialUsed

  end type TMagFieldInp


  !> Magnetic field, including scalar and vector effects
  type TMagField

    !> Is the magnetic field uniform
    logical :: isUniformBField

    !> Are only scalar (F) or also vector (T) field effects included
    logical :: isVectorPotentialUsed

    !> Is the field normal to xy plane
    logical :: isFieldPlaneNormal

    !> z aligned magnetic field
    real(dp) :: Bz

    !> Choice for z-aligned field, Coulomb gage, +/-1 and 0 correspond to Landau gauge with A
    !> aligned x and y or circular (symmetric) gauge respectively
    real(dp) :: beta

    !> Magnetic field in an arbitrary direction, circular/symmetric Coulomb gauge
    real(dp) :: B(3)

    !> Gauge origin
    real(dp) :: origin(3)

  contains

    procedure :: vecPotTransform
    procedure :: vecPotInvTransform1, vecPotInvTransform2
    generic :: vecPotInvTransform => vecPotInvTransform1, vecPotInvTransform2
    procedure :: addShift
    procedure :: isAPresent

  end type TMagField

  !> Magnetic flux quantum (h / e) in atomic units - corresponds to MKS (not CGS Gauss) system
  real(dp), parameter :: Phi0 = 2.0_dp * pi

contains

  !> Initialise magnetic field derived type
  subroutine TMagField_init(this, inp, latVecs, errStatus)

    !> Instance
    type(TMagField), intent(out) :: this

    !> Input data for magnetic fields
    type(TMagFieldInp), intent(in) :: inp

    !> Lattice vectors, if relevant
    real(dp), intent(in) :: latVecs(:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    this%isFieldPlaneNormal = .false.
    this%Bz = 0.0_dp
    this%beta = 0.0_dp
    this%B(:) = 0.0_dp
    this%origin(:) = 0.0_dp

    this%isVectorPotentialUsed = inp%isVectorPotentialUsed

    if (allocated(inp%B)) then
      this%isFieldPlaneNormal = .false.
      this%B(:) = inp%B
      if (this%isVectorPotentialUsed) then
        this%origin(:) = inp%origin
      else
        this%origin(:) = 0.0_dp
      end if
    else if (allocated(inp%nQuanta)) then
      this%isFieldPlaneNormal = .true.
      this%Bz = phi0 / (real(inp%nQuanta, dp) * sqrt(sum(cross3(latVecs(:,1),latVecs(:,2))**2)))
      write(stdOut, *)"Magnetic field ", this%Bz * 1.0E20_dp * hbar * Coulomb__au * AA__Bohr**2, "T"
      this%beta = 0.0_dp ! symmetric gauge for the moment
    else
      @:RAISE_ERROR(errStatus, -1, "Magnetic vector potential incorrectly initialised")
    end if

  end subroutine TMagField_init


  !> Is there effect of the associated vector potential present for this B field
  pure function isAPresent(this)

    !> Magnetic field type container
    class(TMagField), intent(in) :: this

    logical :: isAPresent

    isAPresent = this%isVectorPotentialUsed

  end function isAPresent


  !> Apply a vector potential to transform hamiltonian (and optionally overlap)
  subroutine vecPotTransform(this, ints, coords, img2centCell, orb, iPair, species, nAtom,&
      & nNeighbour, iNeighbour)

    !> Magnetic field type container
    class(TMagField), intent(in) :: this

    !> Integrals
    type(TIntegral), intent(inout) :: ints

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    !> Maping from image to real atoms
    integer,  intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse matrix index
    integer,  intent(in) :: iPair(0:,:)

    !> Species of atoms
    integer,  intent(in) :: species(:)

    !> Number of atoms in central cell
    integer, intent(in) :: nAtom

    !> Number of neighbouring atoms
    integer, intent(in) :: nNeighbour(:)

    !> Neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    real(dp) :: tmpA(3,2), phase
    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
    integer :: iNeigh

    if (.not.this%isVectorPotentialUsed) then
      return
    end if

    @:ASSERT(allocated(ints%hamiltonianGauged) .eqv. allocated(ints%iHamiltonianGauged))
    @:ASSERT(allocated(ints%overlapGauged) .eqv. allocated(ints%iOverlapGauged))

    if (this%isFieldPlaneNormal) then
      if (allocated(ints%iHamiltonian)) then
        do iAt1 = 1, nAtom
          iSp1 = species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          do iNeigh = 0, nNeighbour(iAt1)
            iAt2 = iNeighbour(iNeigh, iAt1)
            iAt2f = img2CentCell(iAt2)
            iSp2 = species(iAt2f)
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh, iAt1)
            tmpA(:,1) = 0.5_dp * (xyPlane(this%Bz, coords(:,iAt1), this%beta)&
                & + xyPlane(this%Bz, coords(:,iAt2), this%beta))
            tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
            phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
            ints%hamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & cos(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)&
                & + sin(phase) * ints%iHamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
            ints%iHamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & cos(phase) * ints%iHamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)&
                & - sin(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
          end do
        end do
      else
        do iAt1 = 1, nAtom
          iSp1 = species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          do iNeigh = 0, nNeighbour(iAt1)
            iAt2 = iNeighbour(iNeigh, iAt1)
            iAt2f = img2CentCell(iAt2)
            iSp2 = species(iAt2f)
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh, iAt1)
            tmpA(:,1) = 0.5_dp * (xyPlane(this%Bz, coords(:,iAt1), this%beta)&
                & + xyPlane(this%Bz, coords(:,iAt2), this%beta))
            tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
            phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
            ints%hamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & cos(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
            ints%iHamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & sin(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
          end do
        end do
      end if
      if (allocated(ints%overlapGauged)) then
        do iAt1 = 1, nAtom
          iSp1 = species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          do iNeigh = 0, nNeighbour(iAt1)
            iAt2 = iNeighbour(iNeigh, iAt1)
            iAt2f = img2CentCell(iAt2)
            iSp2 = species(iAt2f)
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh, iAt1)
            tmpA(:,1) = 0.5_dp * (xyPlane(this%Bz, coords(:,iAt1), this%beta)&
                & + xyPlane(this%Bz, coords(:,iAt2), this%beta))
            tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
            phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
            ints%overlapGauged(iOrig+1:iOrig+nOrb2*nOrb1) =&
                & cos(phase) * ints%overlap(iOrig+1:iOrig+nOrb2*nOrb1)
            ints%iOverlapGauged(iOrig+1:iOrig+nOrb2*nOrb1) =&
                & sin(phase) * ints%overlap(iOrig+1:iOrig+nOrb2*nOrb1)
          end do
        end do
      end if
    else
      if (allocated(ints%iHamiltonian)) then
        do iAt1 = 1, nAtom
          iSp1 = species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          do iNeigh = 0, nNeighbour(iAt1)
            iAt2 = iNeighbour(iNeigh, iAt1)
            iAt2f = img2CentCell(iAt2)
            iSp2 = species(iAt2f)
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh, iAt1)
            tmpA(:,1) = 0.5_dp * (circular(this%B, coords(:,iAt1), this%origin)&
                & + circular(this%B, coords(:,iAt2), this%origin))
            tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
            phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
            ints%hamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & cos(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)&
                & + sin(phase) * ints%iHamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
            ints%iHamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & cos(phase) * ints%iHamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)&
                & - sin(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
          end do
        end do
      else
        do iAt1 = 1, nAtom
          iSp1 = species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          do iNeigh = 0, nNeighbour(iAt1)
            iAt2 = iNeighbour(iNeigh, iAt1)
            iAt2f = img2CentCell(iAt2)
            iSp2 = species(iAt2f)
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh, iAt1)
            tmpA(:,1) = 0.5_dp * (circular(this%B, coords(:,iAt1), this%origin)&
                & + circular(this%B, coords(:,iAt2), this%origin))
            tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
            phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
            ints%hamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & cos(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
            ints%iHamiltonianGauged(iOrig+1:iOrig+nOrb2*nOrb1,:) =&
                & sin(phase) * ints%hamiltonian(iOrig+1:iOrig+nOrb2*nOrb1,:)
          end do
        end do
      end if
      if (allocated(ints%overlapGauged)) then
        do iAt1 = 1, nAtom
          iSp1 = species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          do iNeigh = 0, nNeighbour(iAt1)
            iAt2 = iNeighbour(iNeigh, iAt1)
            iAt2f = img2CentCell(iAt2)
            iSp2 = species(iAt2f)
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh, iAt1)
            tmpA(:,1) = 0.5_dp * (circular(this%B, coords(:,iAt1), this%origin)&
                & + circular(this%B, coords(:,iAt2), this%origin))
            tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
            phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
            ints%overlapGauged(iOrig+1:iOrig+nOrb2*nOrb1) =&
                & cos(phase) * ints%overlap(iOrig+1:iOrig+nOrb2*nOrb1)
            ints%iOverlapGauged(iOrig+1:iOrig+nOrb2*nOrb1) =&
                & sin(phase) * ints%overlap(iOrig+1:iOrig+nOrb2*nOrb1)
          end do
        end do
      end if
    end if

  end subroutine vecPotTransform


  !> Apply inverse of vector potential transform to a sparse matrix
  subroutine vecPotInvTransform1(this, rho, iRho, coords, img2centCell, orb, iPair, species, nAtom,&
      & nNeighbour, iNeighbour)

    !> Magnetic field type container
    class(TMagField), intent(in) :: this

    !> Matrix real part
    real(dp), intent(inout) :: rho(:)

    !> Matrix imaginary part
    real(dp), intent(inout) :: iRho(:)

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    !> Maping from image to real atoms
    integer,  intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse matrix index
    integer,  intent(in) :: iPair(0:,:)

    !> Species of atoms
    integer,  intent(in) :: species(:)

    !> Number of atoms in central cell
    integer, intent(in) :: nAtom

    !> Number of neighbouring atoms
    integer, intent(in) :: nNeighbour(:)

    !> Neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    real(dp) :: tmpA(3,2), phase
    real(dp) :: tmpRe(orb%mOrb**2), tmpIm(orb%mOrb**2)
    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2, iNeigh

    if (.not.this%isVectorPotentialUsed) then
      return
    end if

    if (this%isFieldPlaneNormal) then
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          tmpA(:,1) = 0.5_dp * (xyPlane(this%Bz, coords(:,iAt1), this%beta)&
              & + xyPlane(this%Bz, coords(:,iAt2), this%beta))
          tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
          phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
          tmpRe(:nOrb2*nOrb1) = cos(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1)&
              & - sin(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1)
          tmpIm(:nOrb2*nOrb1) = cos(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1)&
              & + sin(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1)
          rho(iOrig+1:iOrig+nOrb2*nOrb1) = tmpRe(:nOrb2*nOrb1)
          iRho(iOrig+1:iOrig+nOrb2*nOrb1) = tmpIm(:nOrb2*nOrb1)
        end do
      end do
    else
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          tmpA(:,1) = 0.5_dp * (circular(this%B, coords(:,iAt1), this%origin)&
              & + circular(this%B, coords(:,iAt2), this%origin))
          tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
          phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
          tmpRe(:nOrb2*nOrb1) = cos(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1)&
              & - sin(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1)
          tmpIm(:nOrb2*nOrb1) = cos(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1)&
              & + sin(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1)
          rho(iOrig+1:iOrig+nOrb2*nOrb1) = tmpRe(:nOrb2*nOrb1)
          iRho(iOrig+1:iOrig+nOrb2*nOrb1) = tmpIm(:nOrb2*nOrb1)
        end do
      end do
    end if

  end subroutine vecPotInvTransform1


  !> Apply inverse of vector potential transform to a sparse matrix
  subroutine vecPotInvTransform2(this, rho, iRho, coords, img2centCell, orb, iPair, species, nAtom,&
      & nNeighbour, iNeighbour)

    !> Magnetic field type container
    class(TMagField), intent(in) :: this

    !> Matrix real part
    real(dp), intent(inout) :: rho(:,:)

    !> Matrix imaginary part
    real(dp), intent(inout) :: iRho(:,:)

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    !> Maping from image to real atoms
    integer,  intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse matrix index
    integer,  intent(in) :: iPair(0:,:)

    !> Species of atoms
    integer,  intent(in) :: species(:)

    !> Number of atoms in central cell
    integer, intent(in) :: nAtom

    !> Number of neighbouring atoms
    integer, intent(in) :: nNeighbour(:)

    !> Neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    real(dp) :: tmpA(3,2), phase
    real(dp) :: tmpRe(orb%mOrb**2, size(rho, dim=2)), tmpIm(orb%mOrb**2, size(rho, dim=2))
    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2, iNeigh

    if (.not.this%isVectorPotentialUsed) then
      return
    end if

    if (this%isFieldPlaneNormal) then
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          tmpA(:,1) = 0.5_dp * (xyPlane(this%Bz, coords(:,iAt1), this%beta)&
              & + xyPlane(this%Bz, coords(:,iAt2), this%beta))
          tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
          phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
          tmpRe(:nOrb2*nOrb1, :) = cos(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1,:)&
              & - sin(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1,:)
          tmpIm(:nOrb2*nOrb1, :) = cos(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1,:)&
              & + sin(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1,:)
          rho(iOrig+1:iOrig+nOrb2*nOrb1,:) = tmpRe(:nOrb2*nOrb1, :)
          iRho(iOrig+1:iOrig+nOrb2*nOrb1,:) = tmpIm(:nOrb2*nOrb1, :)
        end do
      end do
    else
      do iAt1 = 1, nAtom
        iSp1 = species(iAt1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 0, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh, iAt1)
          tmpA(:,1) = 0.5_dp * (circular(this%B, coords(:,iAt1), this%origin)&
              & + circular(this%B, coords(:,iAt2), this%origin))
          tmpA(:,2) = coords(:,iAt2) - coords(:,iAt1)
          phase = Phi0 * dot_product(tmpA(:,1), tmpA(:,2))
          tmpRe(:nOrb2*nOrb1, :) = cos(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1,:)&
              & - sin(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1,:)
          tmpIm(:nOrb2*nOrb1, :) = cos(phase) * iRho(iOrig+1:iOrig+nOrb2*nOrb1,:)&
              & + sin(phase) * rho(iOrig+1:iOrig+nOrb2*nOrb1,:)
          rho(iOrig+1:iOrig+nOrb2*nOrb1,:) = tmpRe(:nOrb2*nOrb1, :)
          iRho(iOrig+1:iOrig+nOrb2*nOrb1,:) = tmpIm(:nOrb2*nOrb1, :)
        end do
      end do
    end if

  end subroutine vecPotInvTransform2


  !> Construct and add shift potential for scalar potential part of magnetic field
  subroutine addShift(this, shift, iShift, orb, species, coords)

    !> Magnetic field type container
    class(TMagField), intent(in) :: this

    !> block shift from the scalar magnetic potential
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> imaginary block shift from the potential
    real(dp), intent(inout), allocatable :: iShift(:,:,:,:)

    !> Angular momentum information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> list of the species for each atom
    integer, intent(in) :: species(:)

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    integer :: iAt, nAtom, iSpin, nSpin, iSp, iSh, iOrb, nSpecies, iCart, iStart, iEnd, iAng
    complex(dp), allocatable :: Lz(:,:), Lplus(:,:)
    real(dp), allocatable :: SpeciesL(:,:,:,:)
    real(dp):: B(3), A(3),shiftAt

    nAtom = size(shift,dim=3)
    nSpin = size(shift,dim=4)
    nSpecies = maxval(species(1:nAtom))

    B(:) = getB_(this)

    ! spin Zeeman part
    select case(nSpin)
    case(2) ! z aligned electron spins
      do iAt = 1, nAtom
        iSp = species(iAt)
        do iSh = 1, orb%nShell(iSp)
          do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
            shift(iOrb,iOrb,iAt,2) = shift(iOrb,iOrb,iAt,2) - gfac * mu_B * B(3)
          end do
        end do
      end do
    case(4)
      do iSpin = 2, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
              shift(iOrb,iOrb,iAt,iSpin) = shift(iOrb,iOrb,iAt,iSpin) - gfac * mu_B * B(iSpin-1)
            end do
          end do
        end do
      end do
    end select

    if (this%isVectorPotentialUsed) then

      ! Shift contribution from quadratic diamagnetic term
      if (this%isFieldPlaneNormal) then

        do iAt = 1, nAtom
          A(:) = xyPlane(this%Bz, coords(:,iAt), this%beta)
          shiftAt = 0.5_dp * dot_product(A,A)
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
              shift(iOrb,iOrb,iAt,iSpin) = shift(iOrb,iOrb,iAt,iSpin) + shiftAt
            end do
          end do
        end do

      else

        do iAt = 1, nAtom
          A(:) = circular(this%B, coords(:,iAt), this%origin)
          shiftAt = 0.5_dp * dot_product(A,A)
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
              shift(iOrb,iOrb,iAt,iSpin) = shift(iOrb,iOrb,iAt,iSpin) + shiftAt
            end do
          end do
        end do

      end if

    end if

    if (allocated(iShift)) then

      ! Orbital Zeeman part
      allocate(SpeciesL(orb%mOrb,orb%mOrb,3,nSpecies))
      SpeciesL(:,:,:,:) = 0.0_dp
      allocate(Lz(orb%mOrb,orb%mOrb))
      allocate(Lplus(orb%mOrb,orb%mOrb))
      do iSp = 1, nSpecies
        do iSh = 1, orb%nShell(iSp)
          Lz(:,:) = 0.0_dp
          Lplus(:,:) = 0.0_dp
          iAng = orb%angShell(iSh,iSp)
          call getLOperatorsForSpecies(orb, iSp, Lz, Lplus)
          speciesL(orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1, &
              & orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1,1,iSp) &
              & = aimag(Lplus(1:2*iAng+1,1:2*iAng+1))
          speciesL(orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1, &
              & orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1,2,iSp) &
              & = -real(Lplus(1:2*iAng+1,1:2*iAng+1))
          speciesL(orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1, &
              & orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1,3,iSp) &
              & = aimag(Lz(1:2*iAng+1,1:2*iAng+1))
        end do
      end do
      do iAt = 1, nAtom
        iSp = species(iAt)
        do iSh = 1, orb%nShell(iSp)
          iStart = orb%posShell(iSh,iSp)
          iEnd = orb%posShell(iSh+1,iSp)-1
          do iCart = 1, 3
            iShift(iStart:iEnd,iStart:iEnd,iAt,1) = iShift(iStart:iEnd,iStart:iEnd,iAt,1)&
                & - mu_B * B(iCart) * SpeciesL(iStart:iEnd,iStart:iEnd,iCart,iSp)
          end do
        end do
      end do

    end if

  end subroutine addShift


  !> Energy contributions in field
  subroutine getEnergy(this, eAtom, population, orbitalL, orb, species, coords)

    !> Magnetic field type container
    class(TMagField), intent(in) :: this

    real(dp), intent(out) :: eAtom(:)

    !> populations for atoms
    real(dp), intent(in) :: population(:, :, :)

    !> Imaginary part of input Mulliken block charges
    real(dp), intent(in), allocatable :: orbitalL(:, :, :)

    !> Angular momentum information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> list of the species for each atom
    integer, intent(in) :: species(:)

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    integer :: iAt, iSp, nAtom
    real(dp) :: B(3), A(3), shiftAt

    nAtom = size(eAtom)

    eAtom(:) = 0.0_dp

    B(:) = getB_(this)

    ! Spin Zeeman part
    do iAt = 1, nAtom
      iSp = species(iAt)
      eAtom(iAt) = eAtom(iAt) - gfac * mu_B *&
          & dot_product(sum(population(:orb%nOrbSpecies(iSp), iAt, :), dim=1), B)
    end do

    ! Orbital Zeeman part
    if (allocated(orbitalL)) then
      do iAt = 1, nAtom
        iSp = species(iAt)
        eAtom(iAt) = eAtom(iAt) - mu_B *&
            & dot_product(sum(orbitalL(:, :orb%nOrbSpecies(iSp), iAt), dim=2), B)
      end do
    end if

    ! quadratic diamagnetic term
    if (this%isVectorPotentialUsed) then
      if (this%isFieldPlaneNormal) then
         do iAt = 1, nAtom
          A(:) = xyPlane(this%Bz, coords(:,iAt), this%beta)
          shiftAt = 0.5_dp * dot_product(A,A)
          iSp = species(iAt)
          eAtom(iAt) = eAtom(iAt) + shiftAt * sum(population(:orb%nOrbSpecies(iSp), iAt, 1))
        end do
      else
        do iAt = 1, nAtom
          A(:) = circular(this%B, coords(:,iAt), this%origin)
          shiftAt = 0.5_dp * dot_product(A,A)
          iSp = species(iAt)
          eAtom(iAt) = eAtom(iAt) + shiftAt * sum(population(:orb%nOrbSpecies(iSp), iAt, 1))
        end do
      end if
    end if

  end subroutine getEnergy


  ! Internal routines

  !> Returns magnetic field from structure
  pure function getB_(this)

    !> Magnetic field type container
    type(TMagField), intent(in) :: this

    real(dp) :: getB_(3)

    if (this%isFieldPlaneNormal) then
      ! z aligned field
      getB_(:2) = 0.0_dp
      getB_(3) = this%Bz
    else
      getB_(:) = this%B
    end if

  end function getB_

end module dftbp_dftb_magfield
