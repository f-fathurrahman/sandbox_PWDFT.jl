program test
  implicit none
  integer, parameter :: N = 4
  complex(8) :: a(N)
  complex(8) :: b(N)
  complex(8) :: res
  !
  real(8) :: ddot
  complex(8) :: zdotc

  a = cmplx(1.d0, 2.d0, kind=8)
  b = cmplx(1.d0, 2.d0, kind=8)
  !
  write(*,*) 'a = ', a(1)
  write(*,*) 'b = ', b(1)
  res = ddot(2*N, a, 1, b, 1)
  write(*,*) 'ddot(a,b) = ', res
  !
  res = zdotc(N, a, 1, b, 1)
  write(*,*) 'zdotc(a, b) = ', res


  a = cmplx(1.d0, 2.d0, kind=8)
  b = cmplx(1.d0, -2.d0, kind=8)
  !
  write(*,*) 'a = ', a(1)
  write(*,*) 'b = ', b(1)
  res = ddot(2*N, a, 1, b, 1)
  write(*,*) 'ddot(a,b) = ', res
  !
  res = zdotc(N, a, 1, b, 1)
  write(*,*) 'zdotc(a, b) = ', res


end program