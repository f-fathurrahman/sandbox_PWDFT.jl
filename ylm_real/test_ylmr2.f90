PROGRAM test_ylmr2

  IMPLICIT NONE
  INTEGER :: lmax2
  INTEGER :: Ng
  REAL(8), ALLOCATABLE :: G(:,:), GG(:)
  REAL(8), ALLOCATABLE :: ylm(:,:)
  INTEGER :: ig
  INTEGER :: lmaxkb, lm, l, m

  CHARACTER(8) :: arg1, arg2, arg3, arg4
  REAL(8) :: in_g1, in_g2, in_g3

  IF( iargc() /= 4 ) THEN 
    STOP 'Need exactly four arguments Gx,Gy,Gz,lmax'
  ENDIF 

  ! Read arguments as strings
  CALL getarg(1, arg1)
  CALL getarg(2, arg2)
  CALL getarg(3, arg3)
  CALL getarg(4, arg4)

  ! Convert to appropriate variables
  READ(arg1,*) in_g1
  READ(arg2,*) in_g2
  READ(arg3,*) in_g3
  READ(arg4,*) lmaxkb

  Ng = 1 ! Harcoded
  lmax2 = (lmaxkb+1)**2

  WRITE(*,*) 'lmax  = ', lmaxkb
  WRITE(*,*) 'lmax2 = ', lmax2

  ALLOCATE( G(3,Ng) )
  ALLOCATE( GG(Ng) )
  ALLOCATE( ylm(Ng,lmax2) )

  DO ig = 1, Ng
    G(:,ig) = (/ in_g1, in_g2, in_g3 /)
    GG(ig) = G(1,ig)**2 + G(2,ig)**2 + G(3,ig)**2
  ENDDO

  CALL ylmr2( lmax2, Ng, G, GG, ylm )

  lm = 0
  DO l = 0,lmaxkb
    WRITE(*,*)
    lm = lm + 1
    WRITE(*,'(A,I5)', advance='no') 'lm = ', lm
    WRITE(*,'(F18.10)') ylm(1,lm)
    DO m = 1,l
      lm = lm + 1
      WRITE(*,'(A,I5)', advance='no') 'lm = ', lm
      WRITE(*,'(F18.10)') ylm(1,lm)
      lm = lm + 1
      WRITE(*,'(A,I5)', advance='no') 'lm = ', lm
      WRITE(*,'(F18.10)') ylm(1,lm)
    ENDDO 
  ENDDO

  DEALLOCATE( G )
  DEALLOCATE( GG )
  DEALLOCATE( ylm )

END PROGRAM
