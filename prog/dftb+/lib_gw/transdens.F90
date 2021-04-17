!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for the calculation of Mulliken populations representing KS transition
!> densities
module dftbp_TransDens
  use dftbp_accuracy
  use dftbp_blasroutines
  implicit none

  private
  public transQ_init, transQ, transQ_destr

  interface transQ_init
    module procedure transQ_init_real
    ! module procedure TransQ_kpts
  end interface transQ_init

  interface transQ
    module procedure transQ_real
    ! module procedure TransQ_kpts
  end interface transQ

  !> Private module variables (suffixed with "_" for clarity)
  logical :: tTransQinit_ = .false.
  real(dp), allocatable :: sTimEvec_(:,:) ! Overlap times eigenvectors

contains


  !> Initializes transq by evaluating S*C (real case, square multiply)
  subroutine transQ_init_real(SSqrReal, eigenVec)

    !>  Overlap matrix (square format)
    real(dp), intent(in) :: SSqrReal(:,:)

    !>  Eigenvectors of Hamiltonian
    real(dp), intent(in) :: eigenVec(:,:)

    @:ASSERT(.not. tTransQinit_)
    @:ASSERT(all(shape(SSqrReal) == shape(eigenVec)))

    allocate(sTimEvec_(size(SSqrReal, dim=1), size(SSqrReal, dim=2)))

    call symm(sTimEvec_, 'L', SSqrReal, eigenVec, 'L')

    tTransQinit_ = .true.

  end subroutine transQ_init_real


  !> Computes l-shell resolved transition charges between Orbital iOrb and jOrb (real case)
  subroutine transQ_real(nAtom, mAngAtom, iOrb, jOrb, eigenVec, qTrans)

    !> Number of Atoms
    integer, intent(in) :: nAtom

    !> Max. angular momentum per atom
    integer, intent(in) :: mAngAtom(:)

    !> Kohn-Sham orbital index
    integer, intent(in) :: iOrb

    !> Kohn-Sham orbital index
    integer, intent(in) :: jOrb

    !> Eigenvectors of Hamiltonian
    real(dp), intent(in) :: eigenVec(:,:)

    !> Transition charges, l-shell resolved
    real(dp), intent(out) :: qTrans(:)

    integer :: nOrb, iAt, iL, iM, indE, indQ

    @:ASSERT(tTransQinit_)
    @:ASSERT(size(mAngAtom) == nAtom)
    @:ASSERT(size(eigenVec,dim=1) == size(eigenVec,dim=2))

    nOrb=size(eigenVec,dim=1)

    @:ASSERT(iOrb <= nOrb)
    @:ASSERT(jOrb <= nOrb)

    indQ = 0
    indE = 0
    do iAt = 1,nAtom
      do iL = 0, mAngAtom(iAt)
        indQ = indQ + 1
        qTrans(indQ) = 0.0_dp
        do iM = 0, 2 * iL
          indE = indE + 1
          qTrans(indQ) = qTrans(indQ) + 0.5_dp*(&
              &  eigenVec(indE,iOrb)*sTimEvec_(indE,jOrb)&
              & + eigenVec(indE,jOrb)*sTimEvec_(indE,iOrb) )
        enddo
      enddo
    enddo

  end subroutine transQ_real


  !> Deallocates temp array
  subroutine transQ_destr()

    DEALLOCATE(sTimEvec_)

    tTransQinit_ = .false.

  end subroutine transQ_destr

end module dftbp_TransDens
