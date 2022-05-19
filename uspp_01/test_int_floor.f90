program test

  implicit none
  real(8) :: a

  a = 1.1d0
  write(*,*)
  write(*,*) 'a = ', a
  write(*,*) 'int(a) = ', int(a)
  write(*,*) 'floor(a) = ', floor(a)
  write(*,*) 'nint(a) = ', nint(a)

  a = 1.5d0
  write(*,*)
  write(*,*) 'a = ', a
  write(*,*) 'int(a) = ', int(a)
  write(*,*) 'floor(a) = ', floor(a)
  write(*,*) 'nint(a) = ', nint(a)

  a = 1.6d0
  write(*,*)
  write(*,*) 'a = ', a
  write(*,*) 'int(a) = ', int(a)
  write(*,*) 'floor(a) = ', floor(a)
  write(*,*) 'nint(a) = ', nint(a)

end program
