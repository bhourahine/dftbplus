!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains some widely used types (at the moment only TOrbitals)
module commontypes
  use accuracy
  implicit none
  private

  public :: TOrbitals, TDenseMatIndex, TBoundaryConditions


  !> Contains information about the orbitals of the species/atoms in the system
  type TOrbitals

    !> Nr. of shells for each atomic species (nSpecies)
    integer, allocatable :: nShell(:)

    !> Nr. of orbitals for each atomic species (nSpecies)
    integer, allocatable :: nOrbSpecies(:)

    !> Ang. momentum of the a particular l-shell on a particular species (maxval(nShell), nSpecies)
    integer, allocatable :: angShell(:,:)

    !> The shell which contains the given orbital on an atom
    !> (maxval(nOrbSpecies), nSpecies)
    integer, allocatable :: iShellOrb(:,:)

    !> Starting pos. within the atomic block of the each of the shells of each species
    !> (maxval(nShell)+1, nSpecies)
    integer, allocatable :: posShell(:,:)

    !> Max. nr. of shells for any species
    integer :: mShell

    !> Max. nr. of orbitals for any species
    integer :: mOrb

  end type TOrbitals

  !> Contains index information for large dense matrices like the hamiltonian
  type TDenseMatIndex

    !> Total number of orbitals in system.
    integer :: nOrb

    !> Start of orbitals for atoms in dense  H/S matrices
    integer, allocatable :: iDenseStart(:)

  end type TDenseMatIndex

  !> Boundary conditions on the system
  type TBoundaryConditions

    !> cluster geometry
    logical :: tCluster

    !> periodic geometry
    logical :: tPeriodic

    !> Lattice vectors
    real(dp), allocatable :: latVec(:,:)

    !> Lattice vectors at start of calculation, if needed
    real(dp) :: origLatVec(3,3)

  end type TBoundaryConditions

end module commontypes
