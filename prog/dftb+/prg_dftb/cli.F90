!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains command line interface parsing
module CLI
  use accuracy
  use message
  use io
  
  type cliData
    character(lc) :: prefix = ''
    logical :: tPrefix
  end type cliData

  public
  
contains

  subroutine getCLI(arguments)
    type(cliData), intent(out) :: arguments
    integer :: count
    character(len=mc) :: arg

    count = command_argument_count()

    arguments%tPrefix = .false.
    select case (count)
    case(0)
    case(1)      
      call get_command_argument(1, arg)
      arguments%tPrefix = .true.
      arguments%prefix = trim(arg) // '.'
      write(stdOut,"(A)")'Prefix for major files : ' // trim(arguments%prefix)
    case default
      call get_command(arg)
      call error("Excess command line arguments : " // trim(arg))
    end select
  end subroutine getCLI
  
end module CLI
