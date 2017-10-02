!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program dftbplus
  use main, only : runDftbPlus
  use cli
  use inputdata_module, only : inputData
  use formatout, only : printDftbHeader
  use parser, only : parseHsdInput
  use initprogram, only : initProgramVariables
  implicit none

  character(len=*), parameter :: RELEASE_VERSION = '17.1'
  integer, parameter :: RELEASE_YEAR = 2017

  type(cliData) :: arguments
  type(inputData), allocatable :: input

  call printDftbHeader(RELEASE_VERSION, RELEASE_YEAR)
  allocate(input)
  call getCLI(arguments)
  call parseHsdInput(arguments,input)
  call initProgramVariables(arguments,input)
  deallocate(input)
  call runDftbPlus()

end program dftbplus
