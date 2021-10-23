!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for atomic position 2nd derivatives
module dftbp_derivs_hessian
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_constants, only : quaternionName
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_rotatedegen, only : TRotateDegen, TRotateDegen_init
  use dftbp_derivs_staticresponse, only : response
  use dftbp_derivs_staticperturb, only : init_perturbation
  use dftbp_dftb_dftbplusu, only : TDftbU
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : getOnsitePopulation
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : totalShift
  use dftbp_dftb_slakocont, only : TSlakoCont, getCutOff
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_dftbplus_mainio, only : writeDerivBandOut
  use dftbp_dftbplus_initprogram, only : derivVBandOut, userOut
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip, MPI_SUM
#:endif
  use dftbp_extlibs_tblite, only : TTBLite
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_mixer_mixer, only : TMixer
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_parallelks, only : TParallelKS
  implicit none

  private
  public :: q0Derivs

contains

  !> Response with respect to (q=0) displacement of atoms
  subroutine q0Derivs(env, parallelKS, coord, isAutotestWritten, autotestTagFile,&
      & isTagResultsWritten, resultsTagFile, taggedWriter, isBandWritten, fdDetailedOut,&
      & isDetailedWritten, filling, eigvals, tolDegen, eigVecsReal, eigVecsCplx, ham, over, orb,&
      & nAtom, species, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
      & isRespKernelRPA, sccCalc, maxSccIter, sccTol, isSccConvRequired, nMixElements,&
      & nIneqMixElements, iEqOrbitals, tempElec, Ef, tFixEf, spinW, thirdOrd, dftbU, iEqBlockDftbu,&
      & onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC, pChrgMixer, kPoint, kWeight, iCellVec,&
      & cellVec, nEFermi, errStatus, hamiltonianType, skHamCont, selfEgy, skOverCont, tbLite,&
      & isHelical)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Coordinates of all atoms including images
    real(dp), intent(in) :: coord(:,:)

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

    !> Hamiltonian type
    integer, intent(in) :: hamiltonianType

    !> Container for the SlaKo Hamiltonian integrals
    type(TSlakoCont), intent(in), optional :: skHamCont

    !> On-site energies for each species
    real(dp), intent(in), optional :: selfegy(:,:)

    !> Container for the SlaKo overlap integrals
    type(TSlakoCont), intent(in), optional :: skOverCont

    !> Library interface handler
    type(TTBLite), allocatable, optional :: tblite

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    integer :: iS, iK, iAt, jAt, iCart
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

    integer :: nIter, fd
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
          & atomic position perturbation")
    end if
    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    lpAtom: do iAt = 1, nAtom

      write(stdOut,"(1X,A,I0)")'Derivative with respect to position of atom ', iAt

      lpDir: do iCart = 1, 3

        write(stdOut,"(2X,A,A,A)")'Derivative in ', quaternionName(iCart), " direction"

        select case(hamiltonianType)
        case(hamiltonianTypes%dftb)

        case(hamiltonianTypes%xtb)

        end select

        dqOut(:,:,:) = 0.0_dp
        dqIn(:,:,:) = 0.0_dp
        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dqBlockIn(:,:,:,:) = 0.0_dp
          dqBlockOut(:,:,:,:) = 0.0_dp
        end if

        dEi(:,:,:) = 0.0_dp
        call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
            & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc,&
            & sccTol, isSccRequired, nIter, pChrgMixer, nMixElements, nIneqMixElements, dqIn,&
            & dqOut, rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut,&
            & nSpin, maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite,&
            & dqBlockIn, dqBlockOut, eigVals, transform, dEi, dEf, Ef, dHam, idHam, dRho, idRho,&
            & tempElec, tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec, iCellVec,&
            & eigVecsReal, eigVecsCplx, dPsiReal, dPsiCmplx, errStatus, isHelical=isHelical,&
            & coord=coord)
        @:PROPAGATE_ERROR(errStatus)

      #:if WITH_SCALAPACK
        ! Add up and distribute eigenvalue derivatives from each processor
        call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
      #:endif

        if (isAutotestWritten.or.isTagResultsWritten) then
          dEiTmp(:,:,:,iAt) = dEi
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

      end do lpDir

    end do lpAtom

  end subroutine q0Derivs

end module dftbp_derivs_hessian
