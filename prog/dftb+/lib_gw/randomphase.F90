!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to calculate the RPA polarizabilty, dielectric matrix and its inverse,
!> needed in GW calculation.
module dftbp_RandomPhase
  use dftbp_accuracy, only : dp, brdRPA
  use dftbp_blasRoutines
  use dftbp_lapackRoutines
  use dftbp_TransDens
  implicit none

  private
  public buildPolar, buildEpsInv

contains

  !> Computes the real part of the RPA polarizability <P(w)> for purely imaginary w
  subroutine buildPolar(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec, delEigenVal, omegaPPi,&
      & polar)

    !> Number of Atoms
    integer, intent(in) :: nAtom

    !> Number of orbitals
    integer, intent(in) :: nOrb

    !> Size of GW basis
    integer, intent(in) :: nAng

    !> Max. angular momentum per atom
    integer, intent(in) :: mAngAtom(:)

    !> Occupation of MO's (spin restricted)
    real(dp), intent(in) :: filling(:)

    !> Eigenvectors of Hamiltonian
    real(dp), intent(in) :: eigenVec(:,:)

    !> Eigenvalue differences
    real(dp), intent(in) :: delEigenVal(:,:)

    !> Imaginary part of w
    real(dp), intent(in) :: omegaPPi

    !> <P(w)> polarizability, symmetric
    real(dp), intent(out) :: polar(:,:)

    integer :: iOrb, jOrb, iAng, jAng
    real(dp) :: delOcc
    real(dp), allocatable :: qTrans(:)

    @:ASSERT(size(mAngAtom) == nAtom)
    @:ASSERT(size(filling) == nOrb)
    @:ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    @:ASSERT(all(shape(delEigenVal) == (/ nOrb, nOrb /)))
    @:ASSERT(all(shape(polar) == (/ nAng, nAng /)))
    allocate(qTrans(nAng))
    polar(:,:) = 0.0_dp

    do iOrb = 1, nOrb
      do jOrb = iOrb+1, nOrb
        call transQ(nAtom, mAngAtom, iOrb, jOrb, eigenVec, qTrans)
        delOcc = filling(iOrb) - filling(jOrb)
        if(abs(delOcc) >  epsilon(1.0_dp)) then
          do iAng = 1, nAng
            do jAng = iAng, nAng
              polar(jAng,iAng) = polar(jAng,iAng) &
                  &-2.0_dp*delOcc*qTrans(jAng)*qTrans(iAng)*&
                  &(delEigenVal(iOrb,jOrb) + brdRPA)/&
                  &(omegaPPi**2 + (delEigenVal(iOrb,jOrb) + brdRPA)**2)
            enddo
          enddo
        endif
      enddo
    enddo

    do iAng = 1, nAng
      polar(iAng, iAng+1:nAng) =  polar(iAng+1:nAng, iAng)
    enddo

    deallocate(qTrans)

  end subroutine buildPolar


  !> Computes RPA dielectric matrix and its inverse : ([eps(w)]S)^-1 = (1 - [v]<P(w)>)^-1
  subroutine buildEpsInv(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec, delEigenVal, eeGamma,&
      & omegaPPi, epsInv)

    !> Number of Atoms
    integer, intent(in) :: nAtom

    !> Number of orbitals
    integer, intent(in) :: nOrb

    !> Size of GW basis
    integer, intent(in) :: nAng

    !> Max. angular momentum per atom
    integer, intent(in) :: mAngAtom(:)

    !> Occupation of MO's (spin restricted)
    real(dp), intent(in) :: filling(:)

    !> Eigenvectors of Hamiltonian
    real(dp), intent(in) :: eigenVec(:,:)

    !> Eigenvalue differences
    real(dp), intent(in) :: delEigenVal(:,:)

    !> [v]
    real(dp), intent(in) :: eeGamma(:,:)

    !> Imaginary part of w
    real(dp), intent(in) :: omegaPPi

    !> ([eps(w)]S)^-1
    real(dp), intent(out) :: epsInv(:,:)

    integer :: iAng
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: polar(:,:)

    @:ASSERT(size(mAngAtom) == nAtom)
    @:ASSERT(size(filling) == nOrb)
    @:ASSERT(all(shape(eigenVec) == (/ nOrb, nOrb /)))
    @:ASSERT(all(shape(delEigenVal) == (/ nOrb, nOrb /)))
    @:ASSERT(all(shape(eeGamma) == (/ nAng, nAng /)))
    @:ASSERT(all(shape(epsInv) == (/ nAng, nAng /)))

    allocate(polar(nAng,nAng))
    allocate(ipiv(nAng))

    call buildPolar(nAtom, nOrb, nAng, mAngAtom, filling, eigenVec, delEigenVal, omegaPPi, polar)

    epsInv(:,:) = 0.0_dp
    do iAng = 1, nAng
      epsInv(iAng,iAng) = 1.0_dp
    enddo

    ! Compute dielectric matrix: 1 - [v]<P(w)>
    call symm(epsInv, 'L', eeGamma, polar, 'L', -1.0_dp, 1.0_dp)

    ! Take the inverse
    call getrf(epsInv, ipiv)
    call getri(epsInv, ipiv)

    deallocate(ipiv)
    deallocate(polar)

  end subroutine buildEpsInv

end module dftbp_RandomPhase
