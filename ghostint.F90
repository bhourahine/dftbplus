module GWGhostInt
  implicit none
  private

  public :: ghostint

contains
  
  
  subroutine ghostint(intgh,lmax,nType)
    !     Parameter
    real*8     intgh(30,30)
    integer, intent(in) :: nType

    !     Local
    real*8     x(7),e(7)
    integer    i,j,n,m,k

    !     Common block variables
    integer lmax(nType)

    open(34,file='GHOST.DEF',status='old')
    rewind(34)
    read(34,*) (e(m), m = 1,7)  
    read(34,*) (x(m), m = 1,7)  
    close(34)

    !     Order of integrals:
    !     3s4s 3s4p 3p4s 3p4p 3p'4p 3s3d 3p3d

    !     3s4s
    intgh(1,1) = x(1)
    !     3s4p
    intgh(1,2) = x(2)
    intgh(1,3) = x(2)
    intgh(1,4) = x(2)
    !     3p4s
    intgh(2,1) = x(3)
    intgh(3,1) = x(3)
    intgh(4,1) = x(3)
    !     3p4p
    intgh(2,2) = x(4)
    intgh(3,3) = x(4)
    intgh(4,4) = x(4)
    !     3p'4p
    intgh(2,3) = x(5)
    intgh(3,2) = x(5)
    intgh(2,4) = x(5)
    intgh(4,2) = x(5)
    intgh(3,4) = x(5)
    intgh(4,3) = x(5)
    !     3s3d
    intgh(1,5) = x(6)
    intgh(1,6) = x(6)
    intgh(1,7) = x(6)
    intgh(1,8) = x(6)
    intgh(1,9) = x(6)
    !     3p3d
    intgh(2,5) = x(7)
    intgh(3,5) = x(7)
    intgh(4,5) = x(7)
    intgh(2,6) = x(7)
    intgh(3,6) = x(7)
    intgh(4,6) = x(7)
    intgh(2,7) = x(7)
    intgh(3,7) = x(7)
    intgh(4,7) = x(7)
    intgh(2,8) = x(7)
    intgh(3,8) = x(7)
    intgh(4,8) = x(7)
    intgh(2,9) = x(7)
    intgh(3,9) = x(7)
    intgh(4,9) = x(7)

    return
  end subroutine ghostint


end module GWGhostInt
