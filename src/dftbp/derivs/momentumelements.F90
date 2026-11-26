!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for evaluating momentum matrix elements between filled and empty states, see doi:
!! 10.1103/PhysRevB.63.201101 10.1103/PhysRevB.98.115115
module dftbp_derivs_momentumelements
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag, pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_periodic, only : TAuxNeighbourList, TNeighbourList
  use dftbp_dftb_sparse2dense, only : unpackHSdk
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_wrappedintr, only : TWrappedInt2, TWrappedReal1
  use dftbp_math_blasroutines, only : hemm

  implicit none

  private
  public :: momentumElements

contains


  !> Calculate the momentum elements between occupied and virtual states
  subroutine momentumElements(env, pElement, parallelKS, eigvals, eigVecsCplx, ints, nTrans, getIA,&
      & wij, neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, orb, denseDesc,&
      & iSparseStart, img2CentCell, kPoint, kWeight, coord, rCellVecs, cellVec, iCellVec, latVecs,&
      & recVec2p, densityMatrix, hybridXc, dab, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    complex(dp), allocatable, intent(out) :: pElement(:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    integer, intent(in) :: nTrans(:)

    type(TWrappedInt2), intent(in) :: getIA(:)

    type(TWrappedReal1), intent(in) :: wij(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbouring atoms (symmetric version)
    type(TAuxNeighbourList), intent(in), allocatable :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in), allocatable :: nNeighbourCamSym(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Vectors to units cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Unit cell vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2 pi (i.e. the inverse lattice vectors)
    real(dp), intent(in) :: recVec2p(:,:)

    !> Holds real and complex delta density matrices and pointers
    type(TDensityMatrix), intent(in) :: densityMatrix

    !> Hybrid functional, if in use
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> On-site atomic dipole matrix elements, if available
    real(dp), allocatable, intent(in) :: dab(:,:,:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: iAt, iK, iKS, iS, iCart, nOrbs, iTrans, nLocCol, nLocRow, iFil, iEmp, iOrb, ii, jj
    integer :: nAtom, nOrbAt
    complex(dp), allocatable :: work(:,:), work2(:,:), dHamHyb_dk(:,:,:,:), hamHyb(:,:,:)
    complex(dp), allocatable :: gaugeEigVecs(:,:), gaugePhase(:)
    integer, allocatable :: iOrb2Atom(:)
    real(dp), allocatable :: posOrb(:)
    !! K-point in cartesian (absolute) units
    real(dp) :: kCart(3)

  #:if WITH_SCALAPACK

    @:ASSERT(.false.)

  #:endif

    allocate(pElement(3, maxval(nTrans), parallelKS%nLocalKS), source=(0.0_dp,0.0_dp))

    nOrbs = size(eigVecsCplx, dim=1)
    nLocCol = size(eigVecsCplx, dim=1)
    nLocRow = size(eigVecsCplx, dim=2)

    allocate(work(nLocCol, nLocRow))
    allocate(work2, mold=work)
    allocate(gaugeEigVecs, mold=work)
    allocate(gaugePhase(nOrbs))

    nAtom = size(nNeighbourSK)

    !! Map from basis functions to the atom they are centred on
    allocate(iOrb2Atom(nOrbs))
    allocate(posOrb(nOrbs))
    do iAt = 1, nAtom
      ii = denseDesc%iAtomStart(iAt)
      jj = denseDesc%iAtomStart(iAt+1) - 1
      iOrb2Atom(ii:jj) = iAt
    end do

    if (allocated(hybridXc)) then
      ! Derivatives of the non-local exchange contribution with respect to cartesian k
      call hybridXc%getCamHamiltonian_dkpts(env, denseDesc, orb, ints, densityMatrix,&
          & neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, iCellVec, cellVec,&
          & rCellVecs, iSparseStart, img2CentCell, kPoint, kWeight, latVecs, dHamHyb_dk, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      ! The exchange operator is non-local, hence does not commute with the position operator and
      ! contributes to the velocity operator also via i [R, H_cam] on gauge transformation
      call hybridXc%getCamHamiltonian_kpts(env, denseDesc, orb, ints, densityMatrix,&
          & neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, iCellVec, cellVec,&
          & rCellVecs, iSparseStart, img2CentCell, kPoint, kWeight, hamHyb, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    do iKS = 1, parallelKS%nLocalKS

      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      ! Transform into the atom position gauge (convention II) to evaluate optical matrix
      ! elements, i.e. |c~> = V |c> with V = diag(exp(i k . R_atom)). Note that the k-vector is
      ! needed in cartesian (absolute) units here, as it multiplies cartesian atomic positions
      kCart(:) = 2.0_dp * pi * matmul(recVec2p, kPoint(:,iK))
      do iAt = 1, nAtom
        ii = denseDesc%iAtomStart(iAt)
        jj = denseDesc%iAtomStart(iAt+1)-1
        gaugePhase(ii:jj) = exp(imag * dot_product(kCart, coord(:, iAt)))
      end do
      gaugeEigVecs(:,:) = eigVecsCplx(:,:,iKS)
      do iOrb = 1, nOrbs
        gaugeEigVecs(:, iOrb) = gaugePhase * gaugeEigVecs(:, iOrb)
      end do

      do iCart = 1, 3

        ! First term of (7) from 10.1103/PhysRevB.98.115115 since basis is non-orthogonal
        ! Pulay term - e_i dS/dk |psi>
        work(:,:) = cmplx(0,0,dp)
        call unpackHSdk(work, ints%overlap, kCart, neighbourList%iNeighbour, nNeighbourSK,&
            & coord, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)
        call hemm(work2, 'l', work, gaugeEigVecs)
        do iOrb = 1, nOrbs
          work2(:,iOrb) = -eigvals(iOrb, iK, iS) * work2(:,iOrb)
        end do

        ! dH/dk|psi>
        work(:,:) = cmplx(0,0,dp)
        call unpackHSdk(work, ints%hamiltonian(:,iS), kCart, neighbourList%iNeighbour,&
            & nNeighbourSK, coord, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)
        if (allocated(hybridXc)) then
          ! The CAM parts are evaluated in the cell (convention I) gauge, so transform them into
          ! the same gauge as the analytical derivatives above,
          !   dH^at/dk = V (dH^cell/dk + i [R, H^cell]) V^dagger
          ! (only the lower triangle is needed, as the result is used by hemm below)
          posOrb(:) = coord(iCart, iOrb2Atom)
          do iOrb = 1, nOrbs
            work(:, iOrb) = work(:, iOrb) + gaugePhase * conjg(gaugePhase(iOrb))&
                & * (dHamHyb_dk(:, iOrb, iKS, iCart)&
                & + imag * (posOrb(:) - posOrb(iOrb)) * hamHyb(:, iOrb, iKS))
          end do
        end if
        call hemm(work2, 'l', work, gaugeEigVecs, beta=cmplx(1,0,dp))

        do iTrans = 1, nTrans(iKS)
          iFil = getIA(iKS)%data(iTrans, 1)
          iEmp = getIA(iKS)%data(iTrans, 2)
          pElement(iCart, iTrans, iKS) = dot_product(gaugeEigVecs(:,iEmp), work2(:,iFil))
        end do

        if (allocated(dab)) then
          ! On-site dipole elements

          work(:,:) = cmplx(0,0,dp)
          do iAt = 1, nAtom
            ii = denseDesc%iAtomStart(iAt)
            jj = denseDesc%iAtomStart(iAt+1)-1
            nOrbAt = orb%nOrbAtom(iAt)
            work(ii:jj,ii:jj) = dab(:nOrbAt, :nOrbAt, iAt, iCart)
          end do

          call hemm(work2, 'l', work, gaugeEigVecs)

          ! A dipole matrix element d contributes i omega_ai d to the velocity operator. Note, the
          ! overlap weighting of the intra-atomic dipoles (i omega <a|S d|i>) is neglected here,
          ! being of the same order as the approximations in the dipole elements themselves
          do iTrans = 1, nTrans(iKS)
            iFil = getIA(iKS)%data(iTrans, 1)
            iEmp = getIA(iKS)%data(iTrans, 2)
            pElement(iCart, iTrans, iKS) = pElement(iCart, iTrans, iKS)&
                & + imag * wij(iKS)%data(iTrans) * dot_product(gaugeEigVecs(:,iEmp), work2(:,iFil))
          end do

        end if

      end do

    end do

  end subroutine momentumElements


end module dftbp_derivs_momentumelements
