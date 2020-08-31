!> KD-tree implemetation
module nanoflann
  use iso_c_binding
  implicit none
  private

  public :: kdtree_wrap

  ! interface bindings for C++ routines
  interface

    subroutine kdtree(n, radius2, points) bind(C, name='kdtree')
      import :: c_int, c_double
      integer(c_int), intent(in) :: n
      real(c_double), intent(in) :: radius2
      real(c_double), intent(in) :: points(3,n)
    end subroutine kdtree

  end interface

contains

  !> Wrapper  for kd-tree call
  subroutine kdtree_wrap(points, radius)

    !> Point cloud
    real(c_double), intent(in) :: points(:,:)

    !> Search radius
    real(c_double), intent(in) :: radius

    integer(c_int) :: n

    n = size(points, dim=2)
    call kdtree(n, radius**2, points)

  end subroutine kdtree_wrap

end module nanoflann
