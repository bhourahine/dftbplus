!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Constraints on the electronic ground state
module dftbp_elecconstraints
  use dftbp_accuracy, only : dp
  use dftbp_angmomentum, only : getLOperators
  use dftbp_commontypes, only : TOrbitals
  use dftbp_conjgrad, only : TConjGrad, init
  use dftbp_fire, only : TFire, TFire_init
  use dftbp_gdiis, only : TDIIS, init
  use dftbp_geoopt, only : geoOptTypes, TGeoOpt, init, reset, next
  use dftbp_globalenv, only : stdOut
  use dftbp_hsdutils2, only : convertByMul
  use dftbp_hsdutils, only : getChild, getChildValue, detailedError, detailedWarning, getChildren
  use dftbp_hsdutils, only : getSelectedAtomIndices
  use dftbp_lbfgs, only : TLbfgs, TLbfgs_init
  use dftbp_message, only : error
  use dftbp_specieslist, only : readSpeciesList
  use dftbp_steepdesc, only : TSteepDesc, init
  use dftbp_typegeometry, only : TGeometry
  use dftbp_unitconversion, only : chargeUnits
  use dftbp_wrappedintr, only : TWrappedInt1, TWrappedReal1, TWrappedReal2
  use dftbp_xmlf90, only : fnode, fnodeList, string, char, getNodeName, getLength, getItem1
  implicit none

  private
  public :: TElecConstraintInp, TElecConstraint, TElecConstraint_init, elecConstrParser

  type TElecConstraintInp

    !> Group of atoms in a constraint
    type(TWrappedInt1), allocatable :: atomGrp(:)

    !> Constraint targets for atom groups
    real(dp), allocatable :: atomNc(:)

    !> Direction of constraint in (q,m) space
    type(TWrappedReal1), allocatable :: atomSpinDir(:)

    !> Optimisation algorithm
    integer :: iConstrOpt

    !> Derivative tolerance for constraint
    real(dp) :: constrTol

    !> Alpha factor to introduce variation into the space
    real(dp) :: diisAlpha

    !> Generations in the history
    integer :: diisGens

    !> Maximum step for driver
    real(dp) :: maxStep

  end type TElecConstraintInp


  type TElecConstraint

    !> Value of the constraint
    real(dp), allocatable :: Nc(:)

    !> Potential
    real(dp), allocatable :: Vc(:)

    ! Weighting function for constrain

    !> Atom(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAt(:)

    !> Atomic orbital(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAtOrb(:)

    !> Atomic orbital charge/spin quaternion involved in each constrain
    type(TWrappedReal2), allocatable :: wAtSpin(:)

    !> Weighting factor for orbitals in the constraint
    type(TWrappedReal1), allocatable :: w(:)

    !> Optimizer
    type(TGeoOpt), allocatable :: pVcOpt

  contains

    procedure constrain

  end type TElecConstraint

contains

  !> Parse input for constraints
  subroutine elecConstrParser(node, geo, input, isSpinPol, is2Component)

    !> Node to parse
    type(fnode), pointer :: node

    !> Geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> Electronic constraint settings on exit
    type(TElecConstraintInp), intent(out) :: input

    !> Is this spin polarized
    logical, intent(in) :: isSpinPol

    !> Is this a two component (Pauli) calculation
    logical, intent(in) :: is2Component

    type(fnode), pointer :: val, child, child2, child3
    type(fnodeList), pointer :: children
    type(string) :: modifier
    type(string) :: buffer, buffer2
    integer :: iConstr, nConstr

    call getChildValue(node, "Driver", child, "gdiis", child=child2)
    call getNodeName(child, buffer)

    select case(char(buffer))
    case("steepestdescent")
      input%iConstrOpt = geoOptTypes%steepestDesc
    case("conjugategradient")
      input%iConstrOpt = geoOptTypes%conjugateGrad
    case("gdiis")
      input%iConstrOpt = geoOptTypes%diis
      input%diisAlpha = 1.0E-1_dp
      input%diisGens = 8
    case("lbfgs")
      input%iConstrOpt = geoOptTypes%lbfgs
    case("fire")
      input%iConstrOpt = geoOptTypes%fire
      input%maxStep = 1.0_dp
    case default
      call detailedError(child, "Unknown driver '" // char(buffer) // "'")
    end select
    call getChildValue(node, "tolerance", input%constrTol, 1.0E-8_dp)

    call getChildValue(node, "Regions", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read specification for regions of atoms
    call getChildren(child, "Atoms", children)
    nConstr = getLength(children)
    allocate(input%atomGrp(nConstr))
    allocate(input%atomNc(nConstr))
    allocate(input%atomSpinDir(nConstr))
    do iConstr = 1, nConstr
      call getItem1(children, iConstr, child2)
      call getChildValue(child2, "Domain", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species,&
          & input%atomGrp(iConstr)%data)
      call getChildValue(child2, "Nc", input%atomNc(iConstr))
      ! charge only part for the moment
      if (isSpinPol) then
        if (is2Component) then
          allocate(input%atomSpinDir(iConstr)%data(4))
          input%atomSpinDir(iConstr)%data(1) = 1.0_dp
        else
          allocate(input%atomSpinDir(iConstr)%data(2))
          input%atomSpinDir(iConstr)%data(1) = 1.0_dp
        end if
      else
        allocate(input%atomSpinDir(iConstr)%data(1))
        input%atomSpinDir(iConstr)%data(1) = 1.0_dp
      end if
    end do

  end subroutine elecConstrParser


  !> Initialise the constraints
  subroutine TElecConstraint_init(this, input, orb)

    !> Constrain structure instance
    type(TElecConstraint), intent(out) :: this

    !> Input data structure
    type(TElecConstraintInp), intent(in) :: input

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    ! Constrain optimizer algorithms

    !> Steepest descent driver
    type(TSteepDesc), allocatable :: pSteepDesc

    !> Conjugate gradient driver
    type(TConjGrad), allocatable :: pConjGrad

    !> gradient DIIS driver
    type(TDIIS), allocatable :: pDIIS

    !> lBFGS driver for geometry optimisation
    type(TLbfgs), allocatable :: pLbfgs

    !> FIRE driver for geometry optimisation
    type(TFire), allocatable :: pFire

    integer :: iConstr, nConstr, ii, jj, iAt, iOrb, nOrb, nSpin

    nConstr = size(input%atomGrp)

    allocate(this%Vc(nConstr))
    this%Vc = 0.0_dp ! should enable optional initialization of Vc from input

    allocate(this%pVcOpt)
    !input%constrTol
    select case(input%iConstrOpt)
    case(geoOptTypes%steepestDesc)
      call error("Not yet implemented")
    case(geoOptTypes%conjugateGrad)
      allocate(pConjGrad)
      call init(pConjGrad, nConstr, input%constrTol, 5E-1_dp)
      call init(this%pVcOpt, pConjGrad)
    case(geoOptTypes%diis)
      allocate(pDIIS)
      call init(pDIIS, nConstr, input%constrTol, input%diisAlpha, input%diisGens)
      call init(this%pVcOpt, pDIIS)
    case(geoOptTypes%lbfgs)
      allocate(pLbfgs)
      call TLbfgs_init(pLbfgs, nConstr, input%constrTol, 1E-8_dp, 1E-1_dp, 8, .false.,&
          & .false., .true.)
      call init(this%pVcOpt, pLbfgs)
    case (geoOptTypes%fire)
      allocate(pFire)
      call TFire_init(pFire, nConstr, input%constrTol, input%maxStep)
      call init(this%pVcOpt, pFire)
    end select
    call reset(this%pVcOpt, this%Vc)

    this%Nc = input%AtomNc

    allocate(this%wAt(nConstr))
    allocate(this%wAtOrb(nConstr))
    allocate(this%wAtSpin(nConstr))
    allocate(this%w(nConstr))
    do iConstr = 1, nConstr
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do
      allocate(this%wAt(iConstr)%data(nOrb))
      allocate(this%wAtOrb(iConstr)%data(nOrb))
      nSpin = size(input%atomSpinDir(iConstr)%data)
      allocate(this%wAtSpin(iConstr)%data(nOrb,nSpin))
      allocate(this%w(iConstr)%data(nOrb))
      this%w(iConstr)%data(:) = 1.0_dp
      this%wAt(iConstr)%data(:) = 0
      this%wAtOrb(iConstr)%data(:) = 0
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        do jj = 1, orb%nOrbAtom(iAt)
          this%wAt(iConstr)%data(nOrb+jj) = iAt
          this%wAtOrb(iConstr)%data(nOrb+jj) = jj
          this%wAtSpin(iConstr)%data(nOrb+jj,:) = input%atomSpinDir(iConstr)%data
        end do
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do

    end do

  end subroutine TElecConstraint_init


  !> Apply constraint to system
  subroutine constrain(this, shift, qIn, energy, tConverged)

    !> Instance
    class(TElecConstraint), intent(inout) :: this

    !> Potential from constraints
    real(dp), intent(out) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> Energy
    real(dp), intent(in) :: energy

    !> Contribution to free energy functional from constraint(s)
    real(dp) :: deltaW

    !> Derivative of energy functional with respect to Vc
    real(dp), allocatable :: dWdVc(:)

    !> Derivative of energy functional with respect to Nc
    real(dp) :: dWdNc

    !> Gradient convergence achieved
    logical, intent(out) :: tConverged

    integer :: iConstr, nConstr

    shift(:,:,:,:) = 0.0_dp

    nConstr = size(this%wAt)
    allocate(dWdVc(nConstr))
    dWdVc(:) = 0.0_dp
    deltaW = 0.0_dp
    do iConstr = 1, nConstr

      call qConstraint(shift, deltaW, dWdVc(iConstr), this%Vc(iConstr), this%Nc(iConstr),&
          & this%wAt(iConstr)%data, this%wAtOrb(iConstr)%data, this%wAtSpin(iConstr)%data,&
          & this%w(iConstr)%data, qIn)

    end do

    call next(this%pVcOpt, -(energy + deltaW), -dWdVc, this%Vc, tConverged)

    write(stdOut,*)'Delta W',deltaW
    write(stdOut,*)'dW/dVc',dWdVc
    write(stdOut,*)'Vc',this%Vc

  end subroutine constrain


  !> Constraint on number of electrons
  subroutine qConstraint(shift, deltaW, dWdV, Vc, Nc, wAt, wOrb, wSp, w, qIn)

    !> shift to which contribution is appended
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> Free energy contribution from contraint
    real(dp), intent(inout) :: deltaW

    !> Derivative of free energy with respect to potential
    real(dp), intent(out) :: dWdV

    !> Potential / Lagrange multiplier
    real(dp), intent(in) :: Vc

    !> target value
    real(dp), intent(in) :: Nc

    !> atoms in w to be constrained
    integer, intent(in) :: wAt(:)

    !> orbitals in w to be constrained
    integer, intent(in) :: wOrb(:)

    !> spins in w to be constrained
    real(dp), intent(in) :: wSp(:,:)

    !> Weighting factor
    real(dp), intent(in) :: w(:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    integer :: iW, iS, nSp
    real(dp) :: wn

    nSp = size(wSp, dim=2)
    wn = 0.0_dp
    do iS = 1, nSp
      do iW = 1, size(wAt)
        wn = wn + w(iW) * wSp(iW, iS) * qIn(wOrb(iW),wAt(iW),iS)
      end do
    end do

    write(stdOut,*)'wn',wn
    dWdV = wn - Nc
    deltaW = deltaW + Vc * dWdV

    do iS = 1, nSp
      do iW = 1, size(wAt)
        shift(wOrb(iW),wOrb(iW),wAt(iW),iS) = shift(wOrb(iW),wOrb(iW),wAt(iW),iS)&
            & + Vc * w(iW) * wSp(iW, iS)
      end do
    end do

  end subroutine qConstraint


  !> Constraint on an atomic shell spin (non-colinear)
  subroutine constrainSPot(shift, qIn, orb, species, conAt, conSh, Starget, V, vec)

    !> shift to append contribution
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> atom to be constrained
    integer, intent(in) :: conAt

    !> shell of atom to be constrained
    integer, intent(in) :: conSh

    !> target value
    real(dp), intent(in) :: Starget

    !> weight of the constraint
    real(dp), intent(in) :: V

    !> direction of spin
    real(dp), intent(in) :: vec(3)

    integer :: iOrb, nSpin, iSpin
    real(dp) :: Sshell(3), W, vecNorm(3)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)):orb%posShell(conSh+1,species(conAt))-1, &
        & conAt,2:4),dim=1)

    if (sqrt(sum(Sshell**2)) < 1.0E-8_dp) then
      Sshell = Sshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = Sshell  / sqrt(sum(Sshell**2))

    ! Push S towards required value
    w = V * 0.5_dp*(dot_product(Sshell,vecNorm) - Starget)

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)),orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainSPot


  !> Constraint on atomic shell orbital angular momentum
  subroutine constrainLPot(iShift,qBlockSkew, orb, species, conAt, conSh, Ltarget, V, vec)

    !> shift block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of
    !> Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of L
    real(dp), intent(in) :: Ltarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLoperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp
    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    if (sqrt(sum(Lshell**2)) < 1.0E-8_dp) then
      Lshell = Lshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = Lshell / sqrt(sum(Lshell**2))

    ! Push L towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)-Ltarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

  end subroutine constrainLPot


  !> Constraint on atomic shell projection of orbital angular momentum
  subroutine constrainMPot(shift, qIn, iShift, qBlockSkew, orb, species, conAt, conSh, MTarget, V, &
      & vec)

    !> block shift
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> Imaginary block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of m
    real(dp), intent(in) :: MTarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLOperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp

    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do


    ! Push Mj towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)-MTarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

  end subroutine constrainMPot


  !> Constraint on atomic shell total angular momentum
  subroutine constrainJPot(shift, qIn, iShift, qBlockSkew, orb, species, conAt, conSh, Jtarget, V, &
      & vec)
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> Imaginary block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of J
    real(dp), intent(in) :: Jtarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), Sshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLoperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp

    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = -sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          & * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    if ( sqrt(sum((lshell + 0.5_dp*Sshell)**2)) < 1.0E-8_dp) then
      Sshell = Sshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = (lshell + 0.5_dp*Sshell) / sqrt(sum((lshell + 0.5_dp*Sshell)**2))

    ! Push J towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)+ &
        & 0.5_dp*dot_product(Sshell,vecNorm) -Jtarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainJPot


  !> Constraint on atomic shell projection of angular momentum
  subroutine constrainMjPot(shift, qIn, iShift, qBlockSkew, orb, species, conAt, conSh, MjTarget,&
      & V, vec)

    !> block shift
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> charges
    real(dp), intent(in) :: qIn(:,:,:)

    !> Imaginary block shift
    real(dp), intent(inout) :: iShift(:,:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Atom for constraint
    integer, intent(in) :: conAt

    !> Shell for constraint
    integer, intent(in) :: conSh

    !> value of Mj
    real(dp), intent(in) :: MjTarget

    !> strength of constraint
    real(dp), intent(in) :: V

    !> direction of constrain
    real(dp), intent(in) :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), Sshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call getLOperators(iSh, Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,1) = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,2) = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,orb%posShell(conSh,iSp): &
        & orb%posShell(conSh+1,iSp)-1,3) = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp

    ! identity part
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1)
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    ! Push Mj towards required value
    w = V * 0.5_dp*(dot_product(lshell,vecNorm)+ &
        & 0.5_dp*dot_product(Sshell,vecNorm)-MjTarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainMjPot

end module dftbp_elecconstraints
