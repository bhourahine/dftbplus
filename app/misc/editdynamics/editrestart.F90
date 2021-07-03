!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

program editRestart
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut, withMpi
  use dftbp_io_message, only : error
  use dftbp_timedep_dynamicsrestart, only : writeRestartFile, readRestartFile
  use editdynamics_inputdata, only : TInputData
  use editdynamics_parser, only : parseInput
  implicit none

  character(len=*), parameter :: releaseName = ''
  integer, parameter :: releaseYear = 2021

  type(TInputData), allocatable :: input
  type(TEnvironment) :: env

  call initGlobalEnv()
  call TEnvironment_init(env)

  if (withMpi) then
    call error("Electron dynamics edit code does not work with MPI")
  end if

  allocate(input)
  call parseInput(input)
  deallocate(input)

  call destructGlobalEnv()

end program editRestart
