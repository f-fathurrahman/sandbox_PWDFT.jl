program test

implicit none
real(8) :: v(3)
real(8) :: dnrm2

v(:) = (/ 1.d0, 2.d0, 3.d0 /)
write(*,*) 'drnm2 v = ', dnrm2(3, v, 1)

end program