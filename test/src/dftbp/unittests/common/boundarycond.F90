!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("boundarycond")
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_boundarycond, only : boundaryConditions, TBoundaryConditions
  use dftbp_common_boundarycond, only : TBoundaryConditions_init
  use dftbp_common_status, only : TStatus
  implicit none

#:contains

  #:block TEST_FIXTURE("initialise")

    type(TBoundaryConditions) :: bcs
    type(TStatus) :: errStatus
    integer, allocatable :: img2CentCell(:)
    real(dp), allocatable :: coords(:,:)
    real(dp) :: outMatrix(3,3), inMatrix(3,3), matrix(3,3), vec(3), spin(3), theta, a0
    real(dp), parameter :: eye(3,3) = real(reshape([1,0,0,0,1,0,0,0,1],[3,3]), dp)
    integer :: iAt, n, m

  #:contains

    #:block TEST("initfail")
      call TBoundaryConditions_init(bcs, boundaryConditions%cluster, errStatus,&
          & reshape([1.0_dp,0.0_dp,0.0_dp],[3,1]))
      @:ASSERT(errStatus%hasError())
    #:endblock


    #:block TEST("spinspirals")

      m = 9
      n = 9
      a0 = 1.234_dp


      ! spin spiral along [110] direction, rotating pi on translation by [4,4,0]
      theta = pi
      vec = [1,1,0]
      vec(:) = vec / norm2(vec) ! unit vector
      vec(:) = theta * vec / norm2(a0 * [4,4,0])
      call TBoundaryConditions_init(bcs, boundaryConditions%pbc3d, errStatus, reshape(vec,[3,1]))
      @:ASSERT(.not.errStatus%hasError())

      call latticeInit(m, n, 1, a0, img2CentCell, coords)

      ! central cell spin along z
      spin = 1.0_dp * [0,0,1]

      !write(*,*)size(coords, dim=2)
      !write(*,*)
      do iAt = 1, size(coords, dim=2)
        outMatrix(:,:) = bcs%foldOutSpinMatrix(iAt, img2centCell, coords)
        !write(*,"(A,3F12.6,3F12.6)")"H", coords(:,iAt), matmul(outMatrix, spin)

        inMatrix(:,:) = bcs%foldInSpinMatrix(iAt, img2centCell, coords)
        matrix(:,:) = matmul(outMatrix,inMatrix)
        @:ASSERT(all(abs(matrix - eye) < 128_dp * epsilon(0.0_dp)))

      end do

    #:endblock

  #:endblock TEST_FIXTURE


  !> Initialise a cubic supercell
  subroutine latticeInit(l, m, n, a0, img2CentCell, coords)

    integer, intent(in) :: l, m, n
    real(dp), intent(in) :: a0
    integer, allocatable, intent(out) :: img2CentCell(:)
    real(dp), allocatable, intent(out) :: coords(:, :)

    integer :: ii, jj, kk, iAt

    allocate(img2CentCell(l*m*n), source=1)
    allocate(coords(3,l*m*n), source=0.0_dp)

    iAt = 0
    do ii = 0, n-1
      do jj = 0, m-1
        do kk = 0, l-1
          iAt = iAt + 1
          coords(:, iAt) = a0 * [kk,jj,ii]
        end do
      end do
    end do

  end subroutine latticeInit

#:endblock TEST_SUITE


@:TEST_DRIVER()
