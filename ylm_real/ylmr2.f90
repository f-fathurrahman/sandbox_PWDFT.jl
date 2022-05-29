! from QE-6.0

SUBROUTINE ylmr2 (lmax2, ng, g, gg, ylm)
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical
  !     Recipes but avoiding the calculation of factorials that generate
  !     overflow for lmax > 11
  !
  IMPLICIT NONE 
  INTEGER, PARAMETER :: DP = 8
  REAL(DP), PARAMETER :: PI = 4.d0*atan(1.d0)
  REAL(DP), PARAMETER :: FPI = 4.d0*PI
  !
  INTEGER, INTENT(in) :: lmax2, ng
  REAL(DP), INTENT(in) :: g(3, ng), gg(ng)
  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !
  REAL(DP), INTENT(out) :: ylm(ng,lmax2)
  !
  ! local variables
  !
  REAL(DP), PARAMETER :: eps = 1.0d-9
  REAL(DP), ALLOCATABLE :: cost(:), sent(:), phi(:), Q(:,:,:)
  REAL(DP) :: c, gmod
  INTEGER :: lmax, ig, l, m, lm

  IF (ng < 1 .or. lmax2 < 1) RETURN 

  DO lmax = 0, 25
    IF( (lmax+1)**2 == lmax2 ) GOTO 10
  ENDDO
  WRITE(*,*) 'Error ylmr: l > 25 or wrong number of Ylm required, lmax, lmax2 = ', lmax, lmax2
  STOP

10 CONTINUE
  !
  IF (lmax == 0) THEN
    ylm(:,1) =  sqrt (1.d0 / fpi)
    RETURN 
  ENDIF 

  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  allocate(cost(ng), sent(ng), phi(ng), Q(ng,0:lmax,0:lmax) )
  Q(:,:,:) = 0.d0 ! ffr
  !

  DO ig = 1, ng
    gmod = sqrt(gg(ig))
    IF( gmod < eps ) THEN
      cost(ig) = 0.d0
    ELSE 
      cost(ig) = g(3,ig)/gmod
    ENDIF 
    !
    !  beware the arc tan, it is defined modulo pi
    !
    IF( g(1,ig) > eps ) THEN 
      phi(ig) = atan( g(2,ig)/g(1,ig) )
    ELSEIF( g(1,ig) < -eps ) THEN
      phi(ig) = atan( g(2,ig)/g(1,ig) ) + pi
    ELSE 
      phi(ig) = sign( pi/2.d0,g(2,ig) )
    ENDIF 
    sent(ig) = sqrt(max(0d0,1.d0-cost(ig)**2))
     
    WRITE(*,'(1x,A,F18.10)') 'phi  = ', phi(ig)
    WRITE(*,'(1x,A,F18.10)') 'cost = ', cost(ig)
    WRITE(*,'(1x,A,F18.10)') 'sint = ', sent(ig)

  ENDDO
  !
  !  Q(:,l,m) are defined as sqrt((l-m)!/(l+m)!) * P(:,l,m) where
  !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  !
  lm = 0
  
  DO l = 0, lmax
     
    WRITE(*,'(1x,A,I4)', advance='no') 'l = ', l

    c = sqrt (DBLE(2*l+1) / fpi)
    write(*,*) 'c = ', c
    IF( l == 0 ) THEN 
      DO ig = 1, ng
        Q(ig,0,0) = 1.d0
      ENDDO 
    ELSEIF( l == 1 ) THEN 
      DO ig = 1, ng
        Q(ig,1,0) =  cost(ig)
        Q(ig,1,1) = -sent(ig)/sqrt(2.d0)
      ENDDO
    ELSE 
      !
      !  recursion on l for Q(:,l,m)
      !
      DO m = 0, l-2
        DO ig = 1, ng
          Q(ig,l,m) = cost(ig)*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(ig,l-1,m) &
                      - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(ig,l-2,m)
        ENDDO 
      ENDDO 
      !
      DO ig = 1, ng
        Q(ig,l,l-1) = cost(ig) * sqrt(DBLE(2*l-1)) * Q(ig,l-1,l-1)
      ENDDO 
      !
      DO ig = 1, ng
        Q(ig,l,l) = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent(ig)*Q(ig,l-1,l-1)
      ENDDO 
    ENDIF 
    write(*,*) 'sum Q = ', sum(Q(1,:,:))
    
    !
    ! Y_lm, m = 0
    !
    lm = lm + 1
    WRITE(*,*) 'm = 0: lm = ', lm
    DO ig = 1, ng
      ylm(ig, lm) = c * Q(ig,l,0)
    ENDDO 
    !
    DO m = 1, l
      !
      ! Y_lm, m > 0
      !
      lm = lm + 1
      WRITE(*,*) 'm > 0: lm = ', lm, ' m = ', m
      DO ig = 1, ng
        ylm(ig,lm) = c * sqrt(2.d0) * Q(ig,l,m) * cos(m*phi(ig))
      ENDDO 
      !
      ! Y_lm, m < 0
      !
      lm = lm + 1
      WRITE(*,*) 'm < 0: lm = ', lm,  ' m = ', -m
      DO ig = 1, ng
        ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * sin(m*phi(ig))
      ENDDO  
    ENDDO ! m = 1,l
  
  ENDDO 
  !
  DEALLOCATE(cost, sent, phi, Q)
  !
  RETURN
END SUBROUTINE 

