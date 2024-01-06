program test_matmul

implicit none
real(8) :: dtau(3), res(3)
real(8) :: bg(3,3)

dtau(:) = (/ 2.d0, 4.d0, 5.d0 /)
bg(1,:) = (/ 1.0d0, 2.d0, 3.d0 /)
bg(2,:) = (/ 1.0d0, 4.d0, 3.d0 /)
bg(3,:) = (/ 1.1d0, 2.d0, 3.d0 /)


res(:) = matmul(dtau, bg) ! ... expr1
write(*,*) 'res = ', res

res(:) = matmul(bg, dtau)
write(*,*) 'res = ', res

res(:) = matmul(transpose(bg), dtau)  ! ... same as expr1
write(*,*) 'res = ', res

end program