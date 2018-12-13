!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the GW approximation for DFTB
module GWApprox
  use accuracy, only : dp, deltaSelfEnergy
  use globalenv, only : stdOut
  use blasroutines
  use TransDens
  use GammaMat
  use RandomPhase
  use PPModel
  use gwXandXC
  use environment
  implicit none
  private

  public :: gwDriver

contains

  !> Implements the GW approximation for DFTB based on the article
  !> Niehaus et. al., PRA 71, 022508 (2005)
  subroutine gwDriver(env, nAtom, nType, mAngAtom, mAngSpecie, iAtomStart, specie, coord, eigenVec,&
      & eigenVal, SSqrReal, XCSqr, filling , eeU, xTab, eTab, hubbU, qOutput, referenceN0, eHFX,&
      & eXC, eQP)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of Atoms
    integer, intent(in) :: nAtom

    !> Number of different species
    integer, intent(in) :: nType

    !> Max. angular momentum per atom
    integer, intent(in) :: mAngAtom(:)

    !> Max. angular momentum per species
    integer, intent(in) :: mAngSpecie(:)

    !> Indexing array for square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Species list of atomic species for all atom
    integer, intent(in) :: specie(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Eigenvectors of DFTB Hamiltonian
    real(dp), intent(in) :: eigenVec(:,:)

    !> Eigenvalues of DFTB Hamiltonian
    real(dp), intent(in) :: eigenVal(:)

    !> Overlap matrix (square format)
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Matrix of Exchange-Correlation potential
    real(dp), intent(in) :: XCSqr(:,:)

    !> Electron occupation of orbital (spin restr.)
    real(dp), intent(in) :: filling(:)

    !> Electron repulsion integrals (compressed WF)
    real(dp), intent(in) :: eeU(:,:)

    !> Atomic exchange integrals
    real(dp), intent(in) :: xTab(:,:,:)

    !> Electron repulsion integrals (uncompressed WF)
    real(dp), intent(in) :: eTab(:,:,:)

    !> Hubbard values (ee + xc)
    real(dp), intent(in) :: hubbU(:,:)

    !> m-resolved populations
    real(dp), intent(in) :: qOutput(:,:)

    !> l-resolved atomic populations
    real(dp), intent(in) :: referenceN0(:,:)

    !> Hartree-Fock exchange energy per orb.
    real(dp), intent(out) :: eHFX(:)

    !> Exchange-Correlation energy per orb.
    real(dp), intent(out) :: eXC(:)

    !> Quasi particle energies
    real(dp), intent(out) :: eQP(:)

    integer :: nAng, maxM, iOrb, jOrb, nOrb, iPPM
    real(dp), allocatable :: qTrans(:)
    real(dp), allocatable :: delEigenVal(:,:)
    real(dp), allocatable :: eeGamma(:,:)
    real(dp), allocatable :: epsInv(:,:)
    real(dp), allocatable :: wMinV(:,:)
    real(dp), allocatable :: eVecPP(:,:)
    real(dp), allocatable :: ampPP(:)
    real(dp), allocatable :: eVar(:)
    real(dp), parameter :: omegaPPi(2) = [0.0_dp, 0.5_dp]
    complex(dp), allocatable :: frqPP(:)
    complex(dp), allocatable :: eSelfCor(:)
    complex(dp), allocatable :: zRenorm(:)
    complex(dp), allocatable :: delQP(:)

    nAng = sum(mAngAtom) + nAtom
    nOrb = dot_product(mAngAtom+1,mAngAtom+1)

    @:ASSERT(size(mAngAtom) == nAtom)
    @:ASSERT(size(mAngSpecie) == nType)
    @:ASSERT(size(iAtomStart) == nAtom)
    @:ASSERT(size(specie) == nAtom)
    @:ASSERT(all(shape(coord) == (/ 3, nAtom /)))
    @:ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    @:ASSERT(size(eigenVal) == nOrb)
    @:ASSERT(all(shape(SSqrReal) == shape(eigenVec)))
    @:ASSERT(all(shape(XCSqr) == shape(eigenVec)))
    @:ASSERT(size(filling) == nOrb)
    @:ASSERT(all(shape(eeU) == (/ maxval(mAngAtom+1), nType /)))
    @:ASSERT(all(shape(hubbU) == (/ maxval(mAngAtom+1), nType /)))
    maxM = (maxval(mAngAtom)+1)**2
    @:ASSERT(all(shape(eTab) == (/ nType, maxM, maxM /)))
    @:ASSERT(all(shape(xTab) == (/ nType, maxM, maxM /)))
    @:ASSERT(all(shape(qOutput) == (/ maxval(mAngAtom+1)**2, nAtom /)))
    @:ASSERT(all(shape(referenceN0) == (/ maxval(mAngAtom+1), nType /)))
    @:ASSERT(size(eHFX) == nOrb)
    @:ASSERT(size(eXC) == nOrb)
    @:ASSERT(size(eQP) == nOrb)

    allocate(qTrans(nAng))
    call transQ_init(SSqrReal,eigenVec)
    allocate(delEigenVal(nOrb,nOrb))
    allocate(eeGamma(nAng,nAng))
    allocate(epsInv(nAng,nAng))
    allocate(wMinV(nAng,nAng))
    allocate(eVecPP(nAng,nAng))
    allocate(ampPP(nAng))
    allocate(frqPP(nAng))
    allocate(eVar(nOrb))
    allocate(eSelfCor(nOrb))
    allocate(zRenorm(nOrb))
    allocate(delQP(nOrb))

    write(stdOut,*)'       ==> Start of GW calculation'

    do jOrb = 1,nOrb
      do iOrb = 1,jOrb-1
        delEigenVal(iOrb, jOrb) = eigenVal(jOrb) - eigenVal(iOrb)
      enddo
    enddo

    ! Get matrix of the coulomb interaction [v] (equal to gamma)
    ! eeU Hubbard parameter contain no XC contribution
    call gammaM(env, nAtom, mAngAtom, specie, coord, eeU, eeGamma)

    ! Plasmon-Pole-Model with two test frequencies omegaPPi
    do iPPM = 1, 2

      ! ([eps(w)]S)^-1 = (1 - [v]<P(w)>)^-1
      call buildEpsInv(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec, delEigenVal, eeGamma,&
          & omegaPPi(iPPM), epsInv)

      ! wMinV =  epsInv * [v] - [v]
      wMinV(:,:) = eeGamma(:,:)

      call symm(wMinV, 'R', eeGamma, epsInv, 'L', 1.0_dp, -1.0_dp)

      call plasmonPole(nAng, wMinV, omegaPPi(iPPM), eVecPP, ampPP, frqPP, iPPM)

    enddo

    ! PPM parameters plasmon frequencies (frqPP) and amplitudes (ampPP) are now known.
    ! eVecPP are eigenvectors of wMinV

    ! Get Hartree-Fock exchange (eHFX) and v_xc (eXC) energies

    ! Implements Eq. 17 of the PRA paper
    call hfExchange(nOrb, nAng, nAtom, mAngAtom, iAtomStart, &
        & specie, eeGamma, eigenVec, filling, xTab, eHFX)

    ! Implements Eq. 22b
    call dftbXC(env, nOrb, nAng, nAtom, mAngAtom, mAngSpecie, specie, &
        & coord, XCSqr, eigenVec, eTab, hubbU, qOutput, referenceN0, eXC)

    ! Consider energy dependence of self energy through renormalization
    ! Self energy is numerically differentiated (Eq. 3)
    eVar(:) = eigenVal(:) + deltaSelfEnergy

    ! Computation of the correlation part of the self energy (Eq. 21)
    call selfCor(nOrb, nAng, nAtom, mAngAtom, eigenVec, filling, eVecPP, ampPP, frqPP, eigenVal,&
        & eVar, eSelfCor)

    zRenorm(:) = eSelfCor(:)
    eVar(:) = eigenVal(:)
    call selfCor(nOrb, nAng, nAtom, mAngAtom, eigenVec, filling, eVecPP, ampPP, frqPP, eigenVal,&
        & eVar, eSelfCor)

    zRenorm(:) = (zRenorm(:) -  eSelfCor(:))/deltaSelfEnergy
    zRenorm(:) = 1.0_dp/(1.0_dp - zRenorm(:))

    ! Final quasi particle energies (Eq. 2)
    delQP(:) = zRenorm(:)*(eSelfCor(:) + eHFX(:) - eXC(:))
    eQP(:) = eigenVal(:) + real(delQP(:))

    deallocate(delQP)
    deallocate(zRenorm)
    deallocate(eSelfCor)
    deallocate(eVar)
    deallocate(frqPP)
    deallocate(eVecPP)
    deallocate(ampPP)
    deallocate(epsInv)
    deallocate(wMinV)
    deallocate(eeGamma)
    deallocate(delEigenVal)

    call transQ_destr()

    deallocate(qTrans)

    write(stdOut,*)'       ==> GW calculation finished'

  end subroutine gwDriver

end module GWApprox
