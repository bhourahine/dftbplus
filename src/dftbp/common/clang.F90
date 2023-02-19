!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'error.fypp'

!> Inter-operation with C routines
module dftbp_common_clang
  use dftbp_common_status, only : TStatus
  use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_null_char, c_f_pointer, c_associated,&
      & c_size_t
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none

  private
  public :: c_free, c_strlen, c2fStr, fortranChar, handleOutputFileName

  !> Various C routines which need a Fortran binding
  interface

    !> Clean up memory attached to a C pointer
    subroutine c_free(ptr) bind(c, name="free")
      import :: c_ptr
      implicit none
      !> Pointer to nullify
      type(c_ptr), value :: ptr
    end subroutine c_free

    !> strlen for a dynamically allocated string
    integer(c_size_t) function c_strlen(s) bind(C, name='strlen')
      import :: c_size_t, c_ptr
      !> C string
      type(c_ptr), value :: s
    end function c_strlen

  end interface

contains


  !> Process a null terminated C string into a Fortran string, freeing the C string and allocating
  !> the Fortran one. If the string is not associated (NULL) or does not match the (optional)
  !> expected length, throw an error
  subroutine c2fStr(cstr, fstr, status, expectedLen)

    !> Double pointer to string in C
    type(c_ptr), intent(inout) :: cstr

    !> Resulting string in Fortran
    character(:), allocatable, intent(out) :: fstr

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Expected string length, if available
    integer, intent(in), optional :: expectedLen

    character, pointer :: fstring(:)
    integer(c_size_t) :: strlen

    if (.not.c_associated(cstr)) then
      @:RAISE_ERROR(status, -1, "C string is not associated")
    end if
    strlen = c_strlen(cstr)
    if (present(expectedLen)) then
      if (expectedLen /= strlen) then
        @:RAISE_FORMATTED_ERROR(status, -1, "('C string is ', I0,&
            & ' characters long, not the expected ', I0, 'characters')", strlen, expectedLen)
      end if
    end if
    call c_f_pointer(cstr, fstring, [strlen+1])
    allocate(character(strlen) :: fstr)
    fstr = transfer(fstring(:strlen), fstr)
    fstring => null()
    call c_free(cstr)

  end subroutine c2fStr


  !> Converts a 0-char terminated C-type string into a Fortran string.
  function fortranChar(cstring, maxlen)

    !> C-type string as array
    character(kind=c_char), intent(in) :: cstring(*)

    !> Maximal string length. If C-string is longer, it will be chopped.
    integer, intent(in), optional  :: maxlen

    !> Resulting Fortran string
    character(:, kind=c_char), allocatable :: fortranChar

    integer :: ii, maxlen0

    if (present(maxlen)) then
      maxlen0 = maxlen
    else
      maxlen0 = huge(maxlen0) - 1
    end if

    do ii = 1, maxlen0
      if (cstring(ii) == c_null_char) then
        exit
      end if
    end do
    allocate(character(ii - 1) :: fortranChar)
    fortranChar = transfer(cstring(1 : ii - 1), fortranChar)

  end function fortranChar


  !> Handles the optional output file name (which should be a NULL-ptr if not present)
  subroutine handleOutputFileName(outputFileName, outputUnit, tOutputOpened)
    type(c_ptr), intent(in) :: outputFileName
    integer, intent(out) :: outputUnit
    logical, intent(out) :: tOutputOpened

    character(c_char), pointer :: pOutputFileName
    character(:), allocatable :: fortranFileName

    if (c_associated(outputFileName)) then
      call c_f_pointer(outputFileName, pOutputFileName)
      fortranFileName = fortranChar(pOutputFileName)
      open(newunit=outputUnit, file=fortranFileName, action="write")
      tOutputOpened = .true.
    else
      outputUnit = output_unit
      tOutputOpened = .false.
    end if

  end subroutine handleOutputFileName

end module dftbp_common_clang
