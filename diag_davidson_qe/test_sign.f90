program test_sign
  implicit none 
  real(8), parameter :: SMALL=1e-4
  real(8) :: denm

  denm = -1d0
  write(*,'(1x,F18.10)') sign(SMALL, denm)
end