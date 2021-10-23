!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for static linear response derivative calculations using perturbation methods
module dftbp_derivs_staticperturb
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_constants, only : Hartree__eV, quaternionName
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_fermihelper, only : deltamn
  use dftbp_derivs_rotatedegen, only : TRotateDegen, TRotateDegen_init
  use dftbp_derivs_staticresponse, only : response
  use dftbp_dftb_dftbplusu, only : TDftbU
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : getOnsitePopulation
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : totalShift
#:if not WITH_SCALAPACK
  use dftbp_dftb_sparse2dense, only : unpackHS
#:endif
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_dftbplus_mainio, only : writeDerivBandOut
  use dftbp_dftbplus_initprogram, only : derivVBandOut, userOut
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip, MPI_SUM
#:endif
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_mixer_mixer, only : TMixer
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_parallelks, only : TParallelKS
  implicit none

  private
  public :: staticPerturWrtE, polarizabilityKernel, init_perturbation

contains

  !> Static (frequency independent) perturbation at q=0 with respect to an electric field
  subroutine staticPerturWrtE(env, parallelKS, filling, eigvals, tolDegen, eigVecsReal,&
      & eigVecsCplx, ham, over, orb, nAtom, species, neighbourList, nNeighbourSK, denseDesc,&
      & iSparseStart, img2CentCell, coord, sccCalc, maxSccIter, sccTol, isSccConvRequired,&
      & nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef, tFixEf, spinW, thirdOrd, dftbU,&
      & iEqBlockDftbu, onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC, pChrgMixer, kPoint, kWeight,&
      & iCellVec, cellVec, tPeriodic, polarisability, dEi, dqOut, neFermi, dEfdE, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Tolerance for degeneracy between eigenvalues
    real(dp), intent(in) :: tolDegen

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Is this a periodic geometry
    logical, intent(in) :: tPeriodic

    !> Static electric polarisability
    real(dp), intent(out) :: polarisability(:,:)

    !> Derivatives of eigenvalues, if required
    real(dp), allocatable, intent(inout) :: dEi(:,:,:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), allocatable, intent(inout) :: dqOut(:,:,:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Derivative of the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: dEfdE(:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer :: iAt, iCart

    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)
    real(dp), allocatable :: dEiTmp(:,:,:), dEfdETmp(:)

    integer :: ii

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dRho(:,:), idRho(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    logical :: tSccCalc
    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: transform(:)

    if (tPeriodic) then
      @:RAISE_ERROR(errStatus, -1, "Electric field polarizability not currently implemented for&
          & periodic systems")
    end if

    if (tFixEf) then
      @:RAISE_ERROR(errStatus, -1, "Perturbation expressions not currently implemented for fixed&
          & Fermi energy")
    end if

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of electric polarisability'
    write(stdOut,*)

    call init_perturbation(parallelKS, tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill, filling,&
        & ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, transform, rangesep, sSqrReal, over,&
        & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, dRhoOut, dRhoIn,&
        & dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi, eigvals, tempElec, Ef,&
        & kWeight)

    tSccCalc = allocated(sccCalc)

    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    call TPotentials_init(dPotential,orb,nAtom,nSpin,0,0)

    ! If derivatives of eigenvalues are needed
    if (allocated(dEi)) then
      dEi(:,:,:,:) = 0.0_dp
      allocate(dEiTmp(nOrbs, nKPts, nIndepHam))
    end if

    if (any(tMetallic)) then
      allocate(dEfdE(nIndepHam, 3))
      allocate(dEfdETmp(nIndepHam))
    end if

    ! if derivatives of valence wavefunctions needed. Note these will have an arbitrary set of
    ! global phases
    ! if (allocated(eigVecsReal)) then
    !   allocate(dPsiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), nIndepHam, 3))
    ! else
    !   allocate(dPsiCmplx(size(eigvecsCplx,dim=1), size(eigvecsCplx,dim=2), nKpts, nIndepHam, 3))
    ! end if

    dqOut(:,:,:,:) = 0.0_dp

    ! polarisation direction
    ! note, could MPI parallelise over this
    lpCart: do iCart = 1, 3

      dqIn(:,:,:) = 0.0_dp
      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if

      dPotential%extAtom(:,:) = 0.0_dp
      dPotential%extShell(:,:,:) = 0.0_dp
      dPotential%extBlock(:,:,:,:) = 0.0_dp

      ! derivative wrt to electric field as a perturbation
      do iAt = 1, nAtom
        dPotential%extAtom(iAt,1) = coord(iCart,iAt)
      end do
      call totalShift(dPotential%extShell, dPotential%extAtom, orb, species)
      call totalShift(dPotential%extBlock, dPotential%extShell, orb, species)

      if (allocated(dEfdETmp)) then
        dEfdETmp(:) = 0.0_dp
      end if

      dEiTmp(:,:,:) = 0.0_dp

      call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList, nNeighbourSK,&
          & img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc, sccTol,&
          & isSccConvRequired, maxSccIter, pChrgMixer, nMixElements, nIneqMixElements, dqIn,&
          & dqOut(:,:,:,iCart), rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn,&
          & dRhoOut, nSpin, maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite,&
          & dqBlockIn, dqBlockOut, eigVals, transform, dEiTmp, dEfdETmp, Ef, dHam, idHam,&
          & dRho, idRho, tempElec, tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec,&
          & iCellVec, eigVecsReal, eigVecsCplx, dPsiReal, dPsiCmplx, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      if (allocated(dEfdE)) then
        dEfdE(:,iCart) = dEfdETmp
      end if
      if (allocated(dEiTmp)) then
        dEi(:,:,:,iCart) = dEiTmp
      end if

      if (any(tMetallic)) then
        write(stdOut,*)
        write(stdOut,"(A,2E20.12)")'d E_f / d E_'//trim(quaternionName(iCart+1))//':',&
            & dEfdE(:,iCart)
        write(stdOut,*)
      end if

      do ii = 1, 3
        polarisability(ii, iCart) = -sum(sum(dqOut(:,:nAtom,1,iCart),dim=1)*coord(ii,:nAtom))
      end do

    end do lpCart

  #:if WITH_SCALAPACK
    if (allocated(dEi)) then
      call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
    end if
  #:endif

    write(stdOut,*)
    write(stdOut,*)'Static polarisability (a.u.)'
    do iCart = 1, 3
      write(stdOut,"(3E20.12)")polarisability(:, iCart)
    end do
    write(stdOut,*)

  end subroutine staticPerturWrtE


  !> Response with respect to a potential at atomic sites
  subroutine polarizabilityKernel(env, parallelKS, isAutotestWritten, autotestTagFile,&
      & isTagResultsWritten, resultsTagFile, taggedWriter, isBandWritten, fdDetailedOut,&
      & isDetailedWritten, filling, eigvals, tolDegen, eigVecsReal, eigVecsCplx, ham, over, orb,&
      & nAtom, species, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
      & isRespKernelRPA, sccCalc, maxSccIter, sccTol, isSccConvRequired, nMixElements,&
      & nIneqMixElements, iEqOrbitals, tempElec, Ef, tFixEf, spinW, thirdOrd, dftbU, iEqBlockDftbu,&
      & onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC, pChrgMixer, kPoint, kWeight, iCellVec,&
      & cellVec, nEFermi, errStatus, isHelical, coord)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> should regression test data be written
    logical, intent(in) :: isAutotestWritten

    !> Name of output file
    character(*), intent(in) :: autotestTagFile

    !> Produce machine readable results
    logical, intent(in) :: isTagResultsWritten

    !> Name of output file
    character(*), intent(in) :: resultsTagFile

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Should eigenvalue (band) data derivatives be written to disc
    logical, intent(in) :: isBandWritten

    !> File descriptor for the human readable output
    integer, intent(in) :: fdDetailedOut

    !> Is detailed.out being produced
    logical, intent(in) :: isDetailedWritten

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Tolerance for degeneracy between eigenvalues
    real(dp), intent(in) :: tolDegen

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Should the kernel be evaluated at the RPA level (non-SCC) or self-consistent
    logical, intent(in) :: isRespKernelRPA

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    !> Coordinates of all atoms including images
    real(dp), intent(in), optional :: coord(:,:)

    integer :: iS, iK, iAt, jAt
    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    real(dp), allocatable :: dEi(:,:,:), dEiTmp(:,:,:,:)

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dRho(:,:), idRho(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    logical :: tSccCalc
    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    real(dp), allocatable :: dqNetAtom(:), dqNetAtomTmp(:,:)

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: transform(:)

    real(dp), allocatable :: dEf(:), dqOut(:,:,:), dqOutTmp(:,:,:,:)

    integer :: nIter, fd, iRec
    character(mc) :: atLabel

    logical :: isSccRequired

    if (isRespKernelRPA) then
      nIter = 1
      isSccRequired = .false.
    else
      nIter = maxSccIter
      isSccRequired = isSccConvRequired
    end if

    call init_perturbation(parallelKS, tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill, filling,&
        & ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, transform, rangesep, sSqrReal, over,&
        & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, dRhoOut, dRhoIn,&
        & dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi, eigvals, tempElec, Ef,&
        & kWeight)

    if (tFixEf) then
      @:RAISE_ERROR(errStatus, -1, "Perturbation expressions not currently implemented for fixed&
          & Fermi energy")
    end if

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of atomic polarisability kernel'
    write(stdOut,*)

    allocate(dqOut(orb%mOrb, nAtom, nSpin))
    allocate(dqNetAtom(nAtom))
    allocate(dEi(nOrbs, nKpts, nSpin))
    if (isAutotestWritten.or.isTagResultsWritten) then
      allocate(dEiTmp(nOrbs, nKpts, nSpin, nAtom))
      allocate(dqOutTmp(orb%mOrb, nAtom, nSpin, nAtom))
      allocate(dqNetAtomTmp(nAtom, nAtom))
      dEiTmp(:,:,:,:) = 0.0_dp
      dqOutTmp(:,:,:,:) = 0.0_dp
      dqNetAtomTmp(:,:) = 0.0_dp
    end if
    if (any(tMetallic)) then
      allocate(dEf(nIndepHam))
    end if

    tSccCalc = allocated(sccCalc)

    if (allocated(dftbU) .and. allocated(onsMEs)) then
      @:RAISE_ERROR(errStatus, -1, "Onsite corrected and DFTB+U terms currently not compatible for&
          & perturbation")
    end if
    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    ! if derivatives of valence wavefunctions needed. Note these will have an arbitrary set of
    ! global phases
    ! if (allocated(eigVecsReal)) then
    !   allocate(dPsiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), nIndepHam, 3))
    ! else
    !   allocate(dPsiCmplx(size(eigvecsCplx,dim=1), size(eigvecsCplx,dim=2), nKpts, nIndepHam, 3))
    ! end if

    lpAtom: do iAt = 1, nAtom

      write(stdOut,*)'Derivative with respect to potential at atom ', iAt

      if (isAutotestWritten.or.isTagResultsWritten) then
        iRec = iAt
      else
        iRec = 1
      end if

      dqOut(:,:,:) = 0.0_dp
      dqIn(:,:,:) = 0.0_dp
      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if
      dPotential%extAtom(:,:) = 0.0_dp
      dPotential%extAtom(iAt,1) = 1.0_dp

      dPotential%extShell(:,:,:) = 0.0_dp
      dPotential%extBlock(:,:,:,:) = 0.0_dp
      call totalShift(dPotential%extShell, dPotential%extAtom, orb, species)
      call totalShift(dPotential%extBlock, dPotential%extShell, orb, species)

      dEi(:,:,:) = 0.0_dp
      call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList, nNeighbourSK,&
          & img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc, sccTol,&
          & isSccRequired, nIter, pChrgMixer, nMixElements, nIneqMixElements, dqIn, dqOut,&
          & rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut, nSpin,&
          & maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, dqBlockIn,&
          & dqBlockOut, eigVals, transform, dEi, dEf, Ef, dHam, idHam, dRho, idRho, tempElec,&
          & tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec, iCellVec, eigVecsReal,&
          & eigVecsCplx, dPsiReal, dPsiCmplx, errStatus, isHelical=isHelical, coord=coord)
      @:PROPAGATE_ERROR(errStatus)

      #:if WITH_SCALAPACK
        ! Add up and distribute eigenvalue derivatives from each processor
        call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
      #:endif

      if (isAutotestWritten.or.isTagResultsWritten) then
        dEiTmp(:,:,:,iAt) = dEi
      end if

      write(stdOut,*)'Frontier orbital derivatives'
      do iS = 1, nIndepHam
        do iK = 1, nKpts
          write(stdOut,*)dEi(nFilled(iS, iK), iK, iS), dEi(nEmpty(iS, iK), iK, iS)
        end do
      end do


      call getOnsitePopulation(dRho(:,1), orb, iSparseStart, dqNetAtom)
      write(stdOut,*)'Derivatives of Mulliken and on-site (net) populations'
      do jAt = 1, nAtom
        write(stdOut,*)jAt, sum(dqOut(:,jAt,1)), dqNetAtom(jAt)
      end do

      if (isAutotestWritten.or.isTagResultsWritten) then
        dqOutTmp(:,:,:,iAt) = dqOut
        dqNetAtomTmp(:,iAt) = dqNetAtom
      end if

      if (isBandWritten) then
        write(atLabel,"(A,I0)")'ATOM ',iAt
        if (iAt == 1) then
          call writeDerivBandOut(derivVBandOut, dEi, kWeight, preLabel=atLabel)
        else
          call writeDerivBandOut(derivVBandOut, dEi, kWeight, isFileAppended=.true.,&
              & preLabel=atLabel)
        end if
      end if

      if (env%tGlobalLead .and. isDetailedWritten) then
        write(fdDetailedOut, "(A,I0)")'Derivatives wrt. a potential at atom ', iAt
        write(fdDetailedOut,"(1X,A)")'Frontier orbital energy derivatives (a.u.)'
        write(fdDetailedOut,"(1X,A,T14,A,T28,A)")"Spin Kpt","Last filled","First empty"
        do iS = 1, nIndepHam
          do iK = 1, nKpts
            write(fdDetailedOut,"(1X,I2,I4,2F14.6)")iS, iK, dEi(nFilled(iS, iK), iK, iS),&
                & dEi(nEmpty(iS, iK), iK, iS)
          end do
        end do
        write(fdDetailedOut,*)'Atomic population derivatives (a.u.)'
        write(fdDetailedOut, "(1X,A,T10,A,T22,A)")"Atom","Mulliken","On-site"
        do jAt = 1, nAtom
          write(fdDetailedOut, "(I5, 2F12.6)")jAt, sum(dqOut(:,jAt,1)), dqNetAtom(jAt)
        end do
        write(fdDetailedOut,*)
      end if

    end do lpAtom

    if (isAutotestWritten) then
      open(newunit=fd, file=autotestTagFile, action="write", status="old", position="append")
      call taggedWriter%write(fd, tagLabels%dEigenDV, dEiTmp)
      call taggedWriter%write(fd, tagLabels%dqdV, dqOutTmp)
      call taggedWriter%write(fd, tagLabels%dqnetdV, dqNetAtomTmp)
      close(fd)
    end if
    if (isTagResultsWritten) then
      open(newunit=fd, file=resultsTagFile, action="write", status="old", position="append")
      call taggedWriter%write(fd, tagLabels%dEigenDV, dEiTmp)
      call taggedWriter%write(fd, tagLabels%dqdV, dqOutTmp)
      call taggedWriter%write(fd, tagLabels%dqnetdV, dqNetAtomTmp)
      close(fd)
    end if

  end subroutine polarizabilityKernel


  !> Initialise variables for perturbation
  subroutine init_perturbation(parallelKS, tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
      & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, transform, rangesep, sSqrReal,&
      & over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, dRhoOut, dRhoIn,&
      & dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi, eigvals, tempElec, Ef,&
      & kWeight)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Tolerance for degeneracy between eigenvalues
    real(dp), intent(in) :: tolDegen

    !> Number of orbitals
    integer, intent(out) :: nOrbs

    !> Number of k-points
    integer, intent(out) :: nKpts

    !> Number of spin channels
    integer, intent(out) :: nSpin

    !> Number of separate hamiltonians
    integer, intent(out) :: nIndepHam

    !> Maximum occupation of single particle states
    real(dp), intent(out) :: maxFill

    !> Filling of un-perturbed system
    real(dp), intent(in) :: filling(:,:,:)

    !> Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Number of (partly) filled states in each [nIndepHam,kpt]
    integer, intent(out), allocatable :: nFilled(:,:)

    !> First (partly) empty state in each [nIndepHam,kpt]
    integer, intent(out), allocatable :: nEmpty(:,:)

    !> Derivative of hamiltonian
    real(dp), intent(out), allocatable :: dHam(:,:)

    !> Derivative of sparse density matrix
    real(dp), intent(out), allocatable :: dRho(:,:)

    !> Imaginary part of derivative of hamiltonian
    real(dp), intent(out), allocatable :: idHam(:,:)

    !> Imaginary part of derivative of sparse density matrix
    real(dp), intent(out), allocatable :: idRho(:,:)

    !> For orbital transformations in the  case of degeneracies
    type(TRotateDegen), intent(out), allocatable :: transform(:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(in) :: rangeSep

    !> Square matrix for overlap (if needed in range separated calculation)
    real(dp), allocatable, intent(out) :: sSqrReal(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(out) :: dRhoIn(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(out) :: dRhoOut(:)

    !> delta density matrix for rangeseparated calculations
    real(dp), pointer :: dRhoInSqr(:,:,:)

    !> delta density matrix for rangeseparated calculations
    real(dp), pointer :: dRhoOutSqr(:,:,:)

    ! derivative of potentials
    type(TPotentials), intent(out) :: dPotential

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Is there a finite density of states at the Fermi energy
    logical, allocatable, intent(out) :: tMetallic(:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    integer :: iS, iK, iLev, ii

    nOrbs = size(filling,dim=1)
    nKpts = size(filling,dim=2)
    nSpin = size(ham,dim=2)

    select case(nSpin)
    case(1,4)
      nIndepHam = 1
    case(2)
      nIndepHam = 2
    end select
    select case(nSpin)
    case(1)
      maxFill = 2.0_dp
    case(2,4)
      maxFill = 1.0_dp
    end select

    allocate(nFilled(nIndepHam, nKpts))
    allocate(nEmpty(nIndepHam, nKpts))
    nFilled(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKPts
        do iLev = 1, nOrbs
          if ( filling(iLev,iK,iS) < epsilon(1.0) ) then
            ! assumes Fermi filling, so above this is empty
            nFilled(iS,iK) = iLev - 1
            exit
          end if
        end do
        ! check if channel is fully filled
        if (nFilled(iS, iK) < 0) then
          nFilled(iS, iK) = nOrbs
        end if
      end do
    end do
    nEmpty(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKpts
        do iLev = 1, nOrbs
          if ( abs( filling(iLev,iK,iS) - maxFill ) > epsilon(1.0)) then
            ! assumes Fermi filling, so this is filled
            nEmpty(iS, iK) = iLev
            exit
          end if
        end do
        !> Check is channel is empty
        if (nEmpty(iS, iK) < 0) then
          nEmpty(iS, iK) = 1
        end if
      end do
    end do

    allocate(dHam(size(ham,dim=1),nSpin))
    allocate(dRho(size(ham,dim=1),nSpin))
    if (nSpin == 4) then
      allocate(idHam(size(ham,dim=1),nSpin))
      allocate(idRho(size(ham,dim=1),nSpin))
      idHam(:,:) = 0.0_dp
    end if

    allocate(transform(parallelKS%nLocalKS))
    do ii = 1, size(transform)
      call TRotateDegen_init(transform(ii), tolDegen)
    end do

    if (allocated(rangeSep)) then
      allocate(sSqrReal(nOrbs, nOrbs))
      sSqrReal(:,:) = 0.0_dp
    #:if not WITH_SCALAPACK
      call unpackHS(sSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
    #:endif
      allocate(dRhoOut(nOrbs * nOrbs * nSpin))
      dRhoOutSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoOut(:nOrbs*nOrbs*nSpin)
      allocate(dRhoIn(nOrbs * nOrbs * nSpin))
      dRhoInSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoIn(:nOrbs*nOrbs*nSpin)
    else
      dRhoInSqr => null()
      dRhoOutSqr => null()
    end if

    call TPotentials_init(dPotential, orb, nAtom, nSpin, 0,0 )

    allocate(tMetallic(nIndepHam, nKpts))
    tMetallic(:,:) = .not.(nFilled == nEmpty -1)
    if (any(tMetallic)) then
      write(stdOut,*)'Metallic system'
      ! Density of electrons at the Fermi energy, required to correct later for shift in Fermi level
      ! at q=0 in metals
      if (allocated(neFermi)) then
        deallocate(neFermi)
      end if
      allocate(neFermi(nIndepHam))
      do iS = 1, nIndepHam
        neFermi(iS) = 0.0_dp
        do iK = 1, nKpts
          do ii = nEmpty(iS, iK), nFilled(iS, iK)
            neFermi(iS) = neFermi(iS) + kWeight(iK) * deltamn(Ef(iS), eigvals(ii,iK,iS), tempElec)
          end do
        end do
      end do
      neFermi(:) = maxFill * neFermi
      write(stdOut,*)'Density of states at the Fermi energy Nf (a.u.):', neFermi
    else
      write(stdOut,*)'Non-metallic system'
    end if

  end subroutine init_perturbation

end module dftbp_derivs_staticperturb
