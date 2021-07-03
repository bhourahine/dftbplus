!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for parsing input to
module editdynamics_parser
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_unitconversion, only : massUnits
  use dftbp_extlibs_xmlf90, only : fnode, removeChild, string, char, textNodeName, fnodeList,&
      & getLength, getNodeName, getItem1, destroyNodeList, destroyNode, assignment(=)
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils, only : getChild, getChildValue, getChildren, getSelectedAtomIndices,&
      & getSelectedIndices, detailedError, detailedWarning
  use dftbp_io_hsdutils2, only : convertByMul, setUnprocessed, warnUnprocessedNodes, getNodeName2
  use dftbp_io_message, only : error
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_type_linkedlist, only : TListCharLc, TListRealR1, TListString, init, destruct, append,&
      & get, len, asArray
  use dftbp_type_oldskdata, only : TOldSkData, readFromFile
  use dftbp_type_typegeometry, only : TGeometry, reduce, setLattice
  use editdynamics_inputdata, only : TInputData
  implicit none
  private

  !> program version
  character(len=*), parameter :: version =  "0.01"

  !> root node name of the input tree
  character(len=*), parameter :: rootTag = "edit"

  !> input file name
  character(len=*), parameter :: hsdInput = "editdump_in.hsd"

  !> parsed output name
  character(len=*), parameter :: hsdParsedInput = "editdump_pin.hsd"

  !> version of the input document
  integer, parameter :: parserVersion = 1

  public :: parseInput

  !! Variables from the Option block

  !> If the program should be verbose
  logical, public :: isVerbose

contains


  !> Read program settings
  subroutine parseInput(ctrl)

    type(TInputData), intent(out) :: ctrl

    type(TOldSKData) :: skData
    type(fnode), pointer :: root, node, tmp, hsdTree
    type(fnode), pointer :: value, child, child2
    type(TListRealR1) :: realBuffer
    type(string) :: buffer, buffer2
    type(TListString) :: lStr
    integer :: inputVersion
    integer :: ii, iSp1, iAt
    real(dp), allocatable :: speciesMass(:), replacementMasses(:)
    type(TListCharLc), allocatable :: skFiles(:)
    character(lc) :: prefix, suffix, separator, elem1, strTmp, filename
    logical :: tLower, tExist
    logical :: tWriteHSD ! HSD output?

    !> geometry of the system
    type(TGeometry) :: geom

    !! Write header
    write(stdout, "(A)") repeat("=", 80)
    write(stdout, "(A)") "  EDIT EHRENFEST RESTART :" // version
    write(stdout, "(A,/)") repeat("=", 80)

    !! Read in input file as HSD
    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)

    write(stdout, "(A)") "Interpreting input file '" // hsdInput // "'"
    write(stdout, "(A)") repeat("-", 80)

    !! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
          &// i2c(parserVersion) // ") do not match")
    end if

  end subroutine parseInput

end module editdynamics_parser
