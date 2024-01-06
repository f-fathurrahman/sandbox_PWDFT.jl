!-----------------
program test_anint
!-----------------
  
implicit none
real(8) :: x8

x8 = 4.321d0
write(*,*)
write(*,*) 'x8 = ', x8
write(*,*) 'aninit(x8) = ', anint(x8)

x8 = 4.621d0
write(*,*)
write(*,*) 'x8 = ', x8
write(*,*) 'aninit(x8) = ', anint(x8)

end program test_anint