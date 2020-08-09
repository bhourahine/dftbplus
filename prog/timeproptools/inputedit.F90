!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising the time propagation file editor
module dftbp_inputedit
  use dftbp_globalenv, only : stdOut
  use dftbp_accuracy, only : dp, lc
  use dftbp_typegeometryhsd, only : TGeometry, readTGeometryGen, readTGeometryXyz
  use dftbp_hsdparser, only : parseHSD, dumpHSD, getNodeHSDName
  use dftbp_hsdutils, only : getChild, getChildren, getChildValue, detailedError
  use dftbp_xmlf90, only : fnode, fnodelist, string, getLength, getItem1, destroyNodeList, char
  use dftbp_formatout, only : writeGenFormat, writeXYZFormat
  use dftbp_message, only : error
  use dftbp_charmanip, only : unquote
  implicit none

  public :: TInputData, parseHsdInput

  !> container for input data constituents
  type TInputData
    logical :: tInitialized = .false.
    type(TGeometry), allocatable :: geom(:)
    character(lc), allocatable :: fragmentLabel(:)
  end type TInputData

  !> Main HSD input file
  character(len=*), parameter :: hsdInputName = "tp_edit_in.hsd"

  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "tp_edit_in"

  !> Version of the current parser
  integer, parameter :: parserVersion = 8

contains

  !> Parse input from an HSD/XML file
  subroutine parseHsdInput(input)

    !> Returns initialised input variables on exit
    type(TInputData), intent(out) :: input

    type(fnode), pointer :: hsdTree
    type(fnode), pointer :: root, child, child2
    type(fnodeList), pointer :: children
    !type(TListString) :: lStr
    type(string) :: buffer
    integer :: iGeom, nGeoms, nChar

    write(stdOut, "(/, A, /)") "***  Parsing and initializing"

    ! Read in the input
    call parseHSD(rootTag, hsdInputName, hsdTree)
    call getChild(hsdTree, rootTag, root)

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    write(stdout, "(A)") "Interpreting input file '" // hsdInputName // "'"
    write(stdout, "(A)") repeat("-", 80)

    ! Read individual structures
    call getChildren(root, "Structure", children)
    nGeoms = getLength(children)
    if (nGeoms == 0) then
      call detailedError(root,"No structures found in input")
    end if
    allocate(input%fragmentLabel(nGeoms))
    do iGeom = 1, nGeoms
      call getItem1(children, iGeom, child)
      call getChildValue(child, "Name", buffer, child=child2)
      nChar = len(trim(unquote(char(buffer))))
      if (nChar > lc .or. nChar < 1) then
        call detailedError(child2,"File name length")
      end if
      input%fragmentLabel(iGeom) = trim(unquote(char(buffer)))
    end do
    call destroyNodeList(children)
    write(*,*)'Structures',nGeoms
    do iGeom = 1, nGeoms
      write(*,*)iGeom, ':',trim(input%fragmentLabel(iGeom))
    end do

  end subroutine parseHsdInput

end module dftbp_inputedit
