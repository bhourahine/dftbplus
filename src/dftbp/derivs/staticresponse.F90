!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for electronic response calculations
module dftbp_derivs_staticresponse
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_linearresponse, only : dRhoStaticReal, dRhoFermiChangeStaticReal,&
      & dRhoStaticCmplx, dRhoFermiChangeStaticCmplx, dRhoStaticPauli, dRhoFermiChangeStaticPauli
  use dftbp_derivs_rotatedegen, only : TRotateDegen
  use dftbp_dftb_blockpothelper, only : appendBlockReduced
  use dftbp_dftb_dftbplusu, only : TDftbU, TDftbU_init, plusUFunctionals
  use dftbp_dftb_onsitecorrection, only : addOnsShift, onsblock_expand
  use dftbp_dftb_orbitalequiv, only : OrbitalEquiv_reduce, OrbitalEquiv_expand
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : mulliken, denseMulliken, getChargePerShell
  use dftbp_dftb_potentials, only : TPotentials
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : addShift, totalShift
  use dftbp_dftb_spin, only : getSpinShift, ud2qm, qm2ud
  use dftbp_dftb_thirdorder, only : TThirdOrder,  TThirdOrderInp, ThirdOrder_init
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip, MPI_SUM
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : DLEN_, scalafx_getdescriptor
#:endif
  use dftbp_io_message, only : warning
  use dftbp_mixer_mixer, only : TMixer, mix, reset
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_parallelks, only : TParallelKS, TParallelKS_init
  implicit none

  private
  public :: response

contains

  !> Evaluates response, given the external perturbation to the hamiltonian (and overlap)
  subroutine response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
      & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc, sccTol,&
      & isSccConvRequired, maxSccIter, pChrgMixer, nMixElements, nIneqMixElements, dqIn, dqOut,&
      & rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut, nSpin, maxFill,&
      & spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, dqBlockIn, dqBlockOut,&
      & eigVals, transform, dEi, dEfdE, Ef, dHam, idHam,  dRho, idRho, tempElec, tMetallic,&
      & neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec, iCellVec, eigVecsReal, eigVecsCplx,&
      & dPsiReal, dPsiCmplx, errStatus, isHelical, coord, dOver, hamPrime, iHamPrime)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    ! derivative of potentials
    type(TPotentials), intent(inout) :: dPotential

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Derivative of sparse Hamiltonian
    real(dp), intent(inout) :: dHam(:,:)

    !> Derivative of imaginary part of sparse Hamiltonian
    real(dp), intent(inout), allocatable :: idHam(:,:)

    !> Derivative of sparse density matrix
    real(dp), intent(inout) :: dRho(:,:)

    !> Derivative of imaginary part of sparse density matrix
    real(dp), intent(inout), allocatable :: idRho(:,:)

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Is there a finite density of states at the Fermi energy
    logical, intent(in) :: tMetallic(:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> Maximum allowed number of electrons in a single particle state
    real(dp), intent(in) :: maxFill

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable, intent(inout) :: transform(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(inout) :: dRhoIn(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(inout) :: dRhoOut(:)

    !> delta density matrix for rangeseparated calculations
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Levels with at least partial filling
    integer, intent(in) :: nFilled(:,:)

    !> Levels that are at least partially empty
    integer, intent(in) :: nEmpty(:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), intent(inout) :: dqOut(:,:,:)

    real(dp), intent(inout), allocatable :: dEi(:,:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), intent(inout) :: dqIn(:,:,:)

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Derivative of Fermi energy
    real(dp), intent(inout), allocatable :: dEfdE(:)

    !> Derivative of block charges (input)
    real(dp), allocatable, intent(inout) :: dqBlockIn(:,:,:,:)

    !> Derivative of block charges (output)
    real(dp), allocatable, intent(inout) :: dqBlockOut(:,:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Square matrix for overlap (if needed)
    real(dp), allocatable, intent(inout) :: sSqrReal(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Derivative of single particle wavefunctions (real case), if needed
    real(dp), allocatable, intent(inout) :: dPsiReal(:,:,:)

    !> Derivative of single particle wavefunctions (complex case), if needed
    complex(dp), allocatable, intent(inout) :: dPsiCmplx(:,:,:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    !> Coordinates of all atoms including images
    real(dp), intent(in), optional :: coord(:,:)

    !> Derivative of the overlap (if relevant to perturbation)
    real(dp), intent(in), optional :: dOver(:)

    !> Derivative of the hamiltonian for a fixed potential (if relevant to perturbation)
    real(dp), intent(in), optional :: hamPrime(:,:)

    !> Imaginary part of derivative of the hamiltonian for a fixed potential (if relevant to
    !> perturbation)
    real(dp), intent(in), optional :: iHamPrime(:,:)

    logical :: isSccCalc, tConverged
    integer :: iSccIter
    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell, nAtom, nSpin)

    integer :: iAt, iKS, iK, iS, iSh, iSp
    real(dp), allocatable :: dRhoExtra(:,:), idRhoExtra(:,:)
    real(dp) :: dqDiffRed(nMixElements)
    real(dp), allocatable :: imdHam(:,:)

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
  #:endif

    @:ASSERT(present(dOver) .eqv. present(hamPrime))
    @:ASSERT(.not.(.not.present(hamPrime) .and. present(iHamPrime)))
    @:ASSERT(.not.(.not.allocated(idHam) .and. present(iHamPrime)))

    isSccCalc = allocated(sccCalc)

    if (allocated(spinW) .or. allocated(thirdOrd)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if
    if (allocated(thirdOrd)) then
      allocate(atomPot(nAtom, nSpin))
    end if

    if (allocated(idHam) .or. present(iHamPrime)) then
      allocate(imdHam(size(over),nSpin))
    end if

    if (any(tMetallic)) then
      allocate(dRhoExtra(size(over),nSpin))
      if (allocated(imdHam)) then
        allocate(idRhoExtra(size(over),nSpin))
      end if
    end if

    if (isSccCalc) then
        call reset(pChrgMixer, size(dqInpRed))
        dqInpRed(:) = 0.0_dp
        dqPerShell(:,:,:) = 0.0_dp
        if (allocated(rangeSep)) then
          dRhoIn(:) = 0.0_dp
          dRhoOut(:) = 0.0_dp
        end if
      end if

      if (isSccCalc .and. maxSccIter > 1) then
        write(stdOut,"(1X,A,T12,A)")'SCC Iter','Error'
      end if

      lpSCC: do iSccIter = 1, maxSccIter

        dPotential%intAtom(:,:) = 0.0_dp
        dPotential%intShell(:,:,:) = 0.0_dp
        dPotential%intBlock(:,:,:,:) = 0.0_dp

        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dPotential%orbitalBlock(:,:,:,:) = 0.0_dp
        end if

        if (isSccCalc .and. iSCCiter>1) then
          call sccCalc%updateCharges(env, dqIn, orb, species)
          call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
          call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1))
          call sccCalc%getShiftPerL(dPotential%intShell(:,:,1))

          if (allocated(spinW)) then
            call getChargePerShell(dqIn, orb, species, dqPerShell)
            shellPot(:,:,:) = 0.0_dp
            call getSpinShift(shellPot(:,:,2:), dqPerShell(:,:,2:), species, orb, spinW)
            dPotential%intShell(:,:,2:) = dPotential%intShell(:,:,2:) + shellPot(:,:,2:)
          end if

          if (allocated(thirdOrd)) then
            atomPot(:,:) = 0.0_dp
            shellPot(:,:,:) = 0.0_dp
            call thirdOrd%getdShiftdQ(atomPot(:,1), shellPot(:,:,1), species, neighbourList, dqIn,&
                & img2CentCell, orb)
            dPotential%intAtom(:,1) = dPotential%intAtom(:,1) + atomPot(:,1)
            dPotential%intShell(:,:,1) = dPotential%intShell(:,:,1) + shellPot(:,:,1)
          end if

          if (allocated(dftbU)) then
            ! note the derivatives of both FLL and pSIC potentials are the same (i.e. pSIC)
            call dftbU%getDftbUShift(dPotential%orbitalBlock, dqBlockIn, species, orb,&
                & plusUFunctionals%pSIC)
          end if
          if (allocated(onsMEs)) then
            ! onsite corrections
            call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock, dqBlockIn, dummy,&
                & onsMEs, species, orb)
          end if

        end if

        call totalShift(dPotential%intShell,dPotential%intAtom, orb, species)
        call totalShift(dPotential%intBlock,dPotential%intShell, orb, species)
        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%orbitalBlock
        end if
        dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%extBlock

        dHam(:,:) = 0.0_dp
        call addShift(dHam, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
            & iSparseStart, nAtom, img2CentCell, dPotential%intBlock)

        if (present(hamPrime)) then
          dHam(:,:) = dHam + hamPrime
        end if

        if (nSpin > 1) then
          dHam(:,:) = 2.0_dp * dHam
          if (allocated(idHam)) then
            imdHam(:,:) = idHam
          end if
          if (present(iHamPrime)) then
            imdHam(:,:) = imdHam + iHamPrime
          end if
          if (allocated(imdHam)) then
            imdHam(:,:) = 2.0_dp * imdHam
          end if
        end if
        call qm2ud(dHam)
        if (allocated(imdHam)) then
          call qm2ud(imdHam)
        end if

        dRho(:,:) = 0.0_dp
        if (allocated(idRho)) then
          idRho(:,:) = 0.0_dp
        end if

        if (.not. allocated(eigVecsReal) .and. present(dOver)) then
          @:RAISE_ERROR(errStatus, -1, "Currently missing atom position derivatives for this case")
        end if

        ! evaluate derivative of density matrix
        if (allocated(eigVecsReal)) then

          do iKS = 1, parallelKS%nLocalKS

            iS = parallelKS%localKS(2, iKS)

            if (allocated(dRhoOut)) then
              ! replace with case that will get updated in dRhoStaticReal
              dRhoOutSqr(:,:,iS) = dRhoInSqr(:,:,iS)
            end if

            call dRhoStaticReal(env, dHam, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
                & denseDesc, iKS, parallelKS, nFilled(:,1), nEmpty(:,1), eigVecsReal, eigVals, Ef,&
                & tempElec, orb, drho(:,iS), dRhoOutSqr, rangeSep, over, nNeighbourLC,&
                & transform(iKS), species,&
                #:if WITH_SCALAPACK
                & desc,&
                #:endif
                & dEi, dPsiReal, errStatus, isHelical, coord, dOver)
            if (errStatus%hasError()) then
              exit
            end if
          end do
          @:PROPAGATE_ERROR(errStatus)

        elseif (nSpin > 2) then

          do iKS = 1, parallelKS%nLocalKS

            iK = parallelKS%localKS(1, iKS)

            call dRhoStaticPauli(env, dHam, imdHam, neighbourList, nNeighbourSK, iSparseStart,&
                & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigvecsCplx, eigVals, Ef,&
                & tempElec, orb, dRho, idRho, kPoint, kWeight, iCellVec, cellVec, iKS,&
                & transform(iKS), species,&
                #:if WITH_SCALAPACK
                & desc,&
                #:endif
                & dEi, dPsiCmplx, errStatus, isHelical, coord)
            if (errStatus%hasError()) then
              exit
            end if
          end do
          @:PROPAGATE_ERROR(errStatus)

          ! adjustment from Pauli to charge/spin
          dRho(:,:) = 2.0_dp * dRho
          if (allocated(idRho)) then
            idRho(:,:) = 2.0_dp * idRho
          end if

        else

          do iKS = 1, parallelKS%nLocalKS

            iK = parallelKS%localKS(1, iKS)

            call dRhoStaticCmplx(env, dHam, neighbourList, nNeighbourSK, iSparseStart,&
                & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigvecsCplx, eigVals, Ef,&
                & tempElec, orb, dRho, kPoint, kWeight, iCellVec, cellVec, iKS, transform(iKS),&
                & species,&
                #:if WITH_SCALAPACK
                & desc,&
                #:endif
                & dEi, dPsiCmplx, errStatus, isHelical, coord)
            if (errStatus%hasError()) then
              exit
            end if
          end do
          @:PROPAGATE_ERROR(errStatus)

        end if

      #:if WITH_SCALAPACK
        ! Add up and distribute density matrix contributions from each group
        call mpifx_allreduceip(env%mpi%globalComm, dRho, MPI_SUM)
      #:endif

        if (any(tMetallic)) then
          ! correct for Fermi level shift for q=0 fields
          dEfdE(:) = 0.0_dp
          dRhoExtra(:,:) = 0.0_dp
          if (allocated(idRhoExtra)) then
            idRhoExtra(:,:) = 0.0_dp
          end if
          do iKS = 1, parallelKS%nLocalKS
            iK = parallelKS%localKS(1, iKS)
            iS = parallelKS%localKS(2, iKS)

            if (.not.tMetallic(iS,iK)) then
              cycle
            end if

            dqOut(:,:,iS) = 0.0_dp
            call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

            dEfdE(iS) = -sum(dqOut(:, :, iS)) / neFermi(iS)

            if (abs(dEfdE(iS)) > 10.0_dp*epsilon(1.0_dp)) then
              ! Fermi level changes, so need to correct for the change in the number of charges

              if (allocated(eigVecsReal)) then

                ! real case, no k-points
                call dRhoFermiChangeStaticReal(dRhoExtra(:, iS), env, parallelKS, iKS,&
                    & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, dEfdE, Ef,&
                    & nFilled(:,iK), nEmpty(:,iK), eigVecsReal, orb, denseDesc, tempElec, eigVals,&
                    & dRhoOutSqr, species,&
                    #:if WITH_SCALAPACK
                    & desc,&
                    #:endif
                    & isHelical, coord)

              elseif (nSpin > 2) then

                ! two component wavefunction cases
                call dRhoFermiChangeStaticPauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS,&
                    & kPoint, kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK,&
                    & img2CentCell, iSparseStart, dEfdE, Ef, nFilled, nEmpty, eigVecsCplx, orb,&
                    & denseDesc, tempElec, eigVals, species,&
                    #:if WITH_SCALAPACK
                    & desc,&
                    #:endif
                    & errStatus, isHelical, coord)
                @:PROPAGATE_ERROR(errStatus)

              else

                ! Complex case with k-points
                call dRhoFermiChangeStaticCmplx(dRhoExtra, env, parallelKS, iKS, kPoint, kWeight,&
                    & iCellVec, cellVec, neighbourList, nNEighbourSK, img2CentCell, iSparseStart,&
                    & dEfdE, Ef, nFilled, nEmpty, eigVecsCplx, orb, denseDesc, tempElec, eigVals,&
                    & species,&
                    #:if WITH_SCALAPACK
                    & desc,&
                    #:endif
                    & isHelical, coord)

              end if

            end if

          end do

          if (nSpin > 2) then
            ! adjustment from Pauli to charge/spin
            dRhoExtra(:,:) = 2.0_dp * dRhoExtra
            if (allocated(idRhoExtra)) then
              idRhoExtra(:,:) = 2.0_dp * idRhoExtra
            end if
          end if

          #:if WITH_SCALAPACK
          ! Add up and distribute density matrix contribution from each group
          call mpifx_allreduceip(env%mpi%globalComm, dRhoExtra, MPI_SUM)
          #:endif
          dRho(:,:) = dRho + dRhoExtra

        end if

        dRho(:,:) = maxFill * drho
        if (allocated(dRhoOut)) then
          dRhoOut(:) = maxFill * dRhoOut
        end if
        call ud2qm(dRho)

        if (allocated(idRho)) then
          idRho(:,:) = maxFill * drho
          call ud2qm(idRho)
        end if

        dqOut(:,:,:) = 0.0_dp
        do iS = 1, nSpin
          call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
              & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
          if (allocated(dftbU) .or. allocated(onsMEs)) then
            dqBlockOut(:,:,:,iS) = 0.0_dp
            call mulliken(dqBlockOut(:,:,:,iS), over, drho(:,iS), orb, neighbourList%iNeighbour,&
                & nNeighbourSK, img2CentCell, iSparseStart)
          end if
        end do

        if (isSccCalc) then

          if (allocated(rangeSep)) then
            dqDiffRed(:) = dRhoOut - dRhoIn
          else
            dqOutRed(:) = 0.0_dp
            call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, dqOutRed(:nIneqMixElements))
            if (allocated(dftbU)) then
              call AppendBlockReduced(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed)
            end if
            if (allocated(onsMEs)) then
              call AppendBlockReduced(dqBlockOut, iEqBlockOnsite, orb, dqOutRed)
            end if
            dqDiffRed(:) = dqOutRed - dqInpRed
          end if
          sccErrorQ = maxval(abs(dqDiffRed))

          if (maxSccIter > 1) then
            write(stdOut,"(1X,I0,T10,E20.12)")iSCCIter, sccErrorQ
          end if
          tConverged = (sccErrorQ < sccTol)

          if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
            if (iSCCIter == 1) then
              if (allocated(rangeSep)) then
                dRhoIn(:) = dRhoOut
                call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
              else
                dqIn(:,:,:) = dqOut
                dqInpRed(:) = dqOutRed
                if (allocated(dftbU) .or. allocated(onsMEs)) then
                  dqBlockIn(:,:,:,:) = dqBlockOut
                end if
              end if

            else

              if (allocated(rangeSep)) then
                call mix(pChrgMixer, dRhoIn, dqDiffRed)
                call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
              else
                call mix(pChrgMixer, dqInpRed, dqDiffRed)
                #:if WITH_MPI
                ! Synchronise charges in order to avoid mixers that store a history drifting apart
                call mpifx_allreduceip(env%mpi%globalComm, dqInpRed, MPI_SUM)
                dqInpRed(:) = dqInpRed / env%mpi%globalComm%size
                #:endif

                call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, orb, dqIn)

                if (allocated(dftbU) .or. allocated(onsMEs)) then
                  dqBlockIn(:,:,:,:) = 0.0_dp
                  if (allocated(dftbU)) then
                    call dftbU%expandBlock(dqInpRed, iEqBlockDFTBU, orb, dqBlockIn, species,&
                        & orbEquiv=iEqOrbitals)
                  else
                    call Onsblock_expand(dqInpRed, iEqBlockOnSite, orb, dqBlockIn,&
                        & orbEquiv=iEqOrbitals)
                  end if
                end if
              end if

            end if

            if (allocated(rangeSep)) then
              call ud2qm(dqIn)
            end if

          end if
        else
          tConverged = .true.
        end if

        if (tConverged) then
          exit lpSCC
        end if

        if (allocated(spinW)) then
          dqPerShell = 0.0_dp
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              dqPerShell(iSh,iAt,:nSpin) = dqPerShell(iSh,iAt,:nSpin) +&
                  & sum(dqIn(orb%posShell(iSh,iSp): orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
            end do
          end do

        end if

      end do lpSCC

      if (isSccCalc .and. .not.tConverged .and. maxSccIter > 1) then
        if (isSccConvRequired) then
          @:RAISE_ERROR(errStatus, -1, "SCC in perturbation is NOT converged, maximal SCC&
              & iterations exceeded")
        else
          call warning("SCC in perturbation is NOT converged, maximal SCC iterations exceeded")
        end if
      end if

    end subroutine response

end module dftbp_derivs_staticresponse
