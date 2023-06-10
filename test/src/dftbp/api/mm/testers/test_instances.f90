!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Reads in several dftb_in.hsd files and run DFTB+ instances simultaneously
program test_instances
  use dftbplus
  use omp_lib
  use testhelpers, only : writeAutotestTag
  implicit none

  ! double precision real
  integer, parameter :: dp = kind(1.0d0)

  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input

  real(dp) :: merminEnergy
  real(dp) :: gradients(3, 2, 2), grossCharges(2, 2), stressTensor(3, 3, 2)

  !integer :: devNull

  character(:), allocatable :: DftbVersion
  integer :: major, minor, patch

  integer :: ii, instance

  call getDftbPlusBuild(DftbVersion)
  write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
  call getDftbPlusApi(major, minor, patch)
  write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

  ! Note: setting the global standard output to /dev/null to suppress output and run-time error
  ! messages
  ! open(newunit=devNull, file="/dev/null", action="write")

  write(*,*)'THREADS',OMP_GET_MAX_THREADS()

  !$omp parallel do default(shared) private(ii, dftbp, input)
  do instance = 1, OMP_GET_MAX_THREADS()

    !call omp_set_num_threads(1)

    ii = mod(instance, 2)

    print "(A,I0,1X,I0)", 'Instance', instance, ii

    ! call TDftbPlus_init(dftbp, outputUnit=devNull)
    call TDftbPlus_init(dftbp)

    ! You should provide the dftb_in.hsd and skfile found in the test/app/dftb+/non-scc/Si_2/ folder
    select case(ii)
    case(1)
      call dftbp%getInputFromFile("dftb1_in.hsd", input)
    case(2)
      call dftbp%getInputFromFile("dftb2_in.hsd", input)
    end select
    call dftbp%setupCalculator(input)

    print "(A,I0,A,I0)", 'Atoms on instance', instance, ' : ', dftbp%nrOfAtoms()

    call dftbp%getEnergy(merminEnergy)
    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy

    call TDftbPlus_destruct(dftbp)

  end do
  !$omp end parallel do

  ! Write file for internal test system
  !call writeAutotestTag(merminEnergy=merminEnergy(1))

end program test_instances
