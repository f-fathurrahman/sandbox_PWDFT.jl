! Author: Paolo Giannozzi in
! Numerical methods in Quantum Mechanics 
! Corso di Laurea Magistrale in Fisica
!
! Modified by Fadjar Fathurrahman
!---------------------------------------------------------------
PROGRAM hydrogen
!---------------------------------------------------------------
  !
  ! Find a solution with given N, L for an hydrogenic atom
  ! solving the radial Schroedinger equation by Numerov method
  ! Atomic (Ry) units
  !
  IMPLICIT NONE 
  !
  INTEGER, PARAMETER :: dp = selected_real_kind(14,200)
  INTEGER :: Nrmesh
  INTEGER :: n, l, i
  REAL(dp) ::  zeta, zmesh, rmax, xmin, dx, e
  REAL(dp), ALLOCATABLE :: r(:), sqr(:), r2(:), y(:), vpot(:)
  !
  ! initialize atomic charge (Z)
  !
  WRITE(*,'(" Atomic Charge > ")', advance='no')
  READ(*,*) zeta
  !write (*,'(" Atomic charge = ",f10.6)') zeta
  IF( zeta < 1.0_dp) stop 'zeta should be >= 1'
  !
  ! initialize logarithmic mesh
  !
  zmesh = zeta 
  rmax  = 100.0_dp
  xmin  = -8.0_dp
  dx    =  0.01_dp
  !
  Nrmesh = (log(zmesh*rmax) - xmin)/dx
  !
  ALLOCATE( r(0:Nrmesh) )
  ALLOCATE( sqr(0:Nrmesh) )
  ALLOCATE( r2(0:Nrmesh) )
  ALLOCATE( vpot(0:Nrmesh) )
  ALLOCATE( y(0:Nrmesh) )
  !
  CALL do_mesh( zmesh, xmin, dx, Nrmesh, r, sqr, r2 )
  !
  ! initialize the potential
  !
  CALL init_pot( zeta, r, Nrmesh, vpot)
  !
  ! open output file that will contain the wavefunctions
  !
  OPEN(7,file='TEMP_wfc.out',status='unknown',form='formatted')
  
  WRITE(*,'(" n, l > ")', advance='no')
  READ(*,*) n,l
  !n = 3
  !l = 1

  ! Check inputs
  IF( n < 1 ) THEN
    DEALLOCATE(y, vpot, r2, sqr, r)
    STOP
  ELSEIF( n < l+1) THEN
    WRITE(*,*) 'error in main: n < l+1 -> wrong number of nodes '
    DEALLOCATE(y, vpot, r2, sqr, r)
    STOP
  ELSEIF( l < 0) THEN
    WRITE(*,*) 'error in main: l < 0 unphysical '
    DEALLOCATE(y, vpot, r2, sqr, r)
    STOP
  ENDIF

  !
  ! solve the Schroedinger equation in radial coordinates by Numerov method
  !
  CALL solve_sheq(n, l, e, Nrmesh, dx, r, sqr, r2, vpot, zeta, y)
  
  !
  ! write out the eigenvalue energy to be compared with the external potential
  !
  WRITE(*,*) 'Energies in Ry'
  WRITE(*,'(" Numeric =",f15.8,",  Analytic =",f15.8)') e, -(zeta/n)**2
  WRITE(7,"('# ',12x,'r',24x,'R(r)',18x,'X(r)=rR(r)',18x,'Veff')")
  DO i=0,Nrmesh
    WRITE(7,*) r(i),y(i)/sqr(i), y(i)*sqr(i), vpot(i)+l*(l+1)/r2(i)
  ENDDO
  WRITE(7,'(/)')
   
END PROGRAM hydrogen

!
!--------------------------------------------------------------------
SUBROUTINE do_mesh( zmesh, xmin, dx, Nrmesh, r, sqr, r2 )
!--------------------------------------------------------------------
  !
  ! initialize radial grid
  !
  IMPLICIT NONE 
  INTEGER, PARAMETER :: dp = selected_real_kind(14,200)
  INTEGER, INTENT(in) :: Nrmesh
  REAL(dp), INTENT(in) :: zmesh, xmin, dx
  REAL(dp), INTENT(out) :: r(0:Nrmesh), sqr(0:Nrmesh), r2(0:Nrmesh)
  !
  INTEGER :: i
  REAL(dp) :: x
  !
  do i = 0,Nrmesh
    x = xmin + dx * i ! equispaced grid
    r(i) = exp(x)/zmesh ! logirithmic grid
    sqr(i)= sqrt(r(i))
    r2(i) = r(i)*r(i)
  ENDDO
  WRITE(*,'(/" radial grid information:")')
  WRITE(*,'(" dx =",f10.6,", xmin =",f10.6,", zmesh =",f10.6)') dx, xmin, zmesh
  WRITE(*,'(" Nrmesh =",i6,", r(0) =",f10.6,", r(Nrmesh) =",f10.6)') Nrmesh, r(0), r(Nrmesh)
  WRITE(*,*) 
  !
  RETURN
END SUBROUTINE do_mesh


!--------------------------------------------------------------------
SUBROUTINE init_pot( zeta, r, Nrmesh, vpot )
!--------------------------------------------------------------------
  !
  ! initialize potential
  !
  IMPLICIT NONE 
  INTEGER, PARAMETER :: dp = selected_real_kind(14,200)
  INTEGER, INTENT(in) :: Nrmesh
  REAL(dp), INTENT(in) :: zeta, r(0:Nrmesh)
  REAL(dp), INTENT(out):: vpot(0:Nrmesh)
  INTEGER :: i

  OPEN(7,file='TEMP_pot.out',status='unknown',form='formatted')
  WRITE(7,'("#       r             V(r)")')
  DO i = 0,Nrmesh
    vpot(i) = -2.0_dp*zeta/r(i) ! in Ry
    WRITE(7,*) r(i),vpot(i)
  ENDDO
  CLOSE(7)
  
  RETURN
END SUBROUTINE init_pot

!---------------------------------------------------------------------
SUBROUTINE solve_sheq( n, l, e, Nrmesh, dx, r, sqr, r2, vpot, zeta, y )
  !---------------------------------------------------------------------
  !
  ! solve the Schroedinger equation in radial coordinates on a 
  ! logarithmic grid by Numerov method - atomic (Ry) units
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = selected_real_kind(14,200), maxiter=100
  REAL(dp), PARAMETER :: eps=1.0D-10
  INTEGER, INTENT(in) :: Nrmesh, n,l
  REAL(dp), INTENT(in) :: dx, r(0:Nrmesh), sqr(0:Nrmesh), r2(0:Nrmesh), & 
       vpot(0:Nrmesh), zeta
  REAL(dp), INTENT(out) :: e, y(0:Nrmesh)
  
  INTEGER :: i, j, icl, nodes, ncross, kkk
  REAL(dp) :: ddx12, sqlhf, x2l2, ycusp, dfcusp, fac, norm, eup, elw, de
  REAL(dp), ALLOCATABLE :: f(:)
  
  ALLOCATE( f(0:Nrmesh) )
  ddx12 = dx*dx/12.0_dp
  sqlhf = (l + 0.5_dp)**2
  x2l2 = 2*l+2
  !
  ! set (very rough) initial lower and upper bounds to the eigenvalue
  !
  eup = vpot(Nrmesh)
  elw = minval(sqlhf/r2(:) + vpot(:))
  IF( eup-elw < eps ) THEN 
    WRITE(*,*) elw, eup
    WRITE(*,*) 'ERROR: solve_sheq: lower and upper bounds are equal'
    STOP
  ENDIF 
  e = 0.5_dp * (elw + eup)

  DO kkk = 1, maxiter
    !
    ! set up the f-function and determine the position of its last
    ! change of sign
    ! f < 0 (approximately) means classically allowed   region
    ! f > 0         "         "        "      forbidden   "
    !
    icl = -1
    f(0) = ddx12 *( sqlhf + r2(0) * (vpot(0)-e) )
    DO i = 1,Nrmesh
      f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
      !
      ! beware: if f(i) is exactly zero the change of sign is not observed
      ! the following line is a trick to prevent missing a change of sign 
      ! in this unlikely but not impossible case:
      !
      IF( f(i) == 0.0_dp ) f(i) = 1.d-20
      IF( f(i) /= sign(f(i),f(i-1)) ) icl=i
    ENDDO

    IF( icl < 0 .or. icl >= Nrmesh-2 ) THEN
      !
      ! classical turning point not found or too far away
      ! no panic: it may follow from a bad choice of eup in
      ! the first iterations. Update e and eup and re-try
      eup = e
      e = 0.5_dp * (eup + elw)
      CYCLE 
    ENDIF 
    !
    ! f function as required by numerov method
    !
    f = 1.0_dp - f
    y = 0
    !
    ! determination of the wave-function in the first two points 
    ! (asymptotic behaviour - second term depends upon the potential)
    !
    nodes = n - l - 1
    y(0) = r(0)**(l+1) *(1.0_dp - 2.0_dp*zeta*r(0)/x2l2) / sqr(0)
    y(1) = r(1)**(l+1) *(1.0_dp - 2.0_dp*zeta*r(1)/x2l2) / sqr(1)
    !
    ! outward integration, count number of crossings
    !
    ncross = 0
    DO i = 1,icl-1
      y(i+1) = ( (12.0_dp - 10.0_dp*f(i))*y(i) - f(i-1)*y(i-1))/f(i+1)
      IF( y(i) /= sign(y(i),y(i+1)) ) ncross = ncross + 1
    ENDDO
    fac = y(icl) 
    !
    ! check number of crossings
    !
    IF( ncross /= nodes ) THEN 
      IF( ncross > nodes ) THEN 
        eup = e
      ELSE
        elw = e
      ENDIF
      e = 0.5_dp *( eup + elw )
      CYCLE 
    ENDIF 
    !
    ! determination of the wave-function in the last two points 
    ! assuming y(Nrmesh+1) = 0 and y(Nrmesh) = dx
    !
    y(Nrmesh) = dx
    y(Nrmesh-1) = ( 12.0_dp - 10.0_dp*f(Nrmesh) )*y(Nrmesh)/f(Nrmesh-1) 
    !
    ! inward integration 
    !
    DO i = Nrmesh-1,icl+1,-1
      y(i-1) = ( (12.0_dp - 10.0_dp*f(i))*y(i) - f(i+1)*y(i+1) )/f(i-1)
      IF( y(i-1) > 1.0d10 ) THEN 
        DO j = Nrmesh,i-1,-1
          y(j) = y(j)/y(i-1)
        ENDDO 
      ENDIF 
    ENDDO 
    !
    ! rescale function to match at the classical turning point (icl)
    !
    fac = fac/y(icl)
    y(icl:) = y(icl:)*fac
    !
    ! normalization - note the change of variable:
    !  \int f(r)dr => \sum_i f_i r_i Delta x
    !
    norm = 0.d0
    DO i=0,Nrmesh
      norm = norm + y(i)*y(i) * r2(i) * dx
    ENDDO 
    norm = sqrt(norm)
    y = y / norm
    
    !
    ! find the value of the cusp at the matching point (icl)
    !
    i = icl
    ycusp = ( y(i-1)*f(i-1) + f(i+1)*y(i+1) + 10.0_dp*f(i)*y(i)) / 12.0_dp
    dfcusp = f(i)*( y(i)/ycusp - 1.0_dp )
    
    !
    ! eigenvalue update using perturbation theory
    !
    de = dfcusp/ddx12 * ycusp*ycusp * dx 
    IF( de > 0.0_dp ) elw = e
    IF( de < 0.0_dp ) eup = e
    
    !
    ! prevent e to go out of bounds, i.e. e > eup or e < elw 
    ! (might happen far from convergence)
    !
    e = max( min (e+de,eup),elw)
    
    !
    ! convergence check
    !
    write(*,'(1x,A,ES18.10)') 'de = ', de
    IF( abs(de) < eps ) EXIT 
    !
  ENDDO 
  !
  ! was convergence achieved ?
  !
  IF( abs(de) > 1.d-10 ) then
    IF( ncross /= nodes ) then
      WRITE(*,*) e, elw, eup, ncross, nodes, icl
    ELSE
      WRITE(*,*) e, de
    ENDIF
    WRITE(*,*) ' error in solve_sheq: too many iterations'
    STOP 
  ELSE
    ! ---- convergence has been achieved -----
    WRITE(*,'(" convergence achieved at iter #",i3," de = ", e10.4)') kkk,de
  ENDIF
  RETURN
END SUBROUTINE solve_sheq


