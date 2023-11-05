!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block MPI_TEST_SUITE("matrixoperations")
  use dftbp_common_accuracy, only : dp
  use dftbp_extlibs_scalapackfx, only : blacsgrid
  use dftbp_math_matrixoperations, only : triangleCopySquareMatrix, blockTriangleCopyHS
  use mpi
  implicit none

  integer :: mycomm, myrank

#:contains

  #:block MPI_TEST_SUITE_INIT
    integer :: error

    mycomm = MPI_COMM_WORLD
    call mpi_comm_rank(mycomm, myrank, error)
    @:MPI_ASSERT(error == 0)

  #:endblock

  #:block MPI_TEST("blockSymmetryMPI")

     integer :: ii, jj, kk

     @:MPI_ASSERT(.false.)

  #:endblock MPI_TEST

#:endblock MPI_TEST_SUITE

#:block MPI_TEST_DRIVER()
#:endblock MPI_TEST_DRIVER
