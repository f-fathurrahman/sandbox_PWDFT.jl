!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )

!SUBROUTINE my_cegterg_v61( &
!  npw, npwx, nvec, nvecx, npol, evc, ethr, &
!  uspp, e, btype, notcnv, lrot, dav_iter )

!
! npol arguments are removed to be able to use matmul
!----------------------------------------------------------------------------
SUBROUTINE my_cegterg_v61( &
  H, S, & 
  npw, npwx, nvec, nvecx, evc, ethr, &
  uspp, e, btype, notcnv, lrot, dav_iter )
!----------------------------------------------------------------------------

!
! ... iterative solution of the eigenvalue problem:
!
! ... ( H - e S ) * evc = 0
!
! ... where H is an hermitean operator, e is a real scalar,
! ... S is an overlap matrix, evc is a complex vector
!

  IMPLICIT NONE
  
  integer, parameter :: DP=8
  
  ! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: npw

  complex(8), intent(in) :: H(npwx,npwx)
  complex(8), intent(in) :: S(npwx,npwx)

  ! leading dimension of matrix evc, as declared in the calling pgm unit
  integer, intent(in) :: npwx

  ! integer number of searched low-lying roots
  integer, intent(in) :: nvec

  ! maximum dimension of the reduced basis set :
  !    (the basis set is refreshed when its dimension would exceed nvecx)
  integer, intent(in) :: nvecx
  
  ! number of spin polarizations
  !integer, intent(in) :: npol

  ! energy threshold for convergence :
  !   root improvement is stopped, when two consecutive estimates of the root
  !   differ by less than ethr.  
  REAL(DP), INTENT(IN) :: ethr

  ! if .FALSE. : do not calculate S|psi>
  LOGICAL, INTENT(IN) :: uspp

  ! band type ( 1 = occupied, 0 = empty )
  INTEGER, INTENT(IN) :: btype(nvec)

  ! .TRUE. if the wfc have already been rotated  
  LOGICAL, INTENT(IN) :: lrot

  ! evc contains the  refined estimates of the eigenvectors  
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)

  ! contains the estimated roots.
  REAL(DP), INTENT(OUT) :: e(nvec)

  ! integer number of iterations performed  
  INTEGER, INTENT(OUT) :: dav_iter
  
  ! number of unconverged roots
  integer, intent(out) :: notcnv


  
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 2
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters

  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  !
  REAL(DP), EXTERNAL :: ddot
  integer :: i

  !
  if( nvec > nvecx / 2 ) then
    write(*,*) 'Error in my_cegterg: nvecx is too small, nvecx = ', nvecx
    write(*,*) 'nvec = ', nvec
    stop
  endif

  ! threshold for empty bands
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  

  ! npol is fixed to 1
  !IF( npol == 1 ) THEN
     kdim = npw
     kdmx = npwx
  !ELSE
  !   kdim = npwx*npol
  !   kdmx = npwx*npol
  !END IF

  ALLOCATE( psi(npwx, nvecx) )
  ALLOCATE( hpsi(npwx, nvecx) )
  IF( uspp ) THEN
    ALLOCATE( spsi(npwx, nvecx) )
  ENDIF

  !
  ALLOCATE( sc(nvecx, nvecx) )
  ALLOCATE( hc(nvecx, nvecx) )  
  ALLOCATE( vc(nvecx, nvecx) )
  ALLOCATE( ew(nvecx) )
  ALLOCATE( conv(nvec) )

  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.

  write(*,*) 'notcnv = ', notcnv

  IF( uspp ) spsi = ZERO

  hpsi = ZERO
  psi  = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  
  ! hpsi contains h times the basis vectors
  ! spsi contains s times the basis vectors
  !CALL h_psi( npwx, npw, nvec, psi, hpsi )
  !IF( uspp ) CALL s_psi( npwx, npw, nvec, psi, spsi )
  
  hpsi(:,1:nvec) = matmul( H, psi(:,1:nvec) )
  if(uspp) spsi(:,1:nvec) = matmul( S, psi(:,1:nvec) )
  
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  
  !
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
              psi, kdmx, hpsi, kdmx, ZERO, hc, nvecx )

  !
  IF( uspp ) THEN
    CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                psi, kdmx, spsi, kdmx, ZERO, sc, nvecx )
  ELSE
    CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                psi, kdmx, psi, kdmx, ZERO, sc, nvecx )
  ENDIF
  

  write(*,*) 'lrot = ', lrot
  IF( lrot ) THEN
    DO n = 1, nbase
      e(n) = REAL( hc(n,n) )
      vc(n,n) = ONE
    END DO
  ELSE
    ! diagonalize the reduced hamiltonian
    CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
    e(1:nvec) = ew(1:nvec)
  ENDIF

  write(*,*) 'After first cdiaghg:'
  do i = 1,nvec
    write(*,'(1x,I4,F18.10)') i, e(i)
  enddo

  !stop 'ffr 207'
  
  ! ----------------------
  ! iterate
  ! ----------------------
  !
  iterate: DO kter = 1, maxter
    !
    write(*,*) 'kter = ', kter
    dav_iter = kter
    !
    np = 0
    !
    DO n = 1, nvec
      IF( .NOT. conv(n) ) THEN
        ! this root not yet converged ... 
        np = np + 1
        ! reorder eigenvectors so that coefficients for unconverged
        ! roots come first. This allows to use quick matrix-matrix 
        ! multiplications to set a new basis vector (see below)
        IF( np /= n ) vc(:,np) = vc(:,n)
        ! for use in g_psi
        ew(nbase+np) = e(n)
      ENDIF
    ENDDO
    !

    nb1 = nbase + 1
    write(*,*) 'nbase = ', nbase
    write(*,*) 'nb1 = ', nb1

    ! expand the basis set with new basis vectors ( H - e*S )|psi> ...
    IF ( uspp ) THEN
       CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, spsi, &
                   kdmx, vc, nvecx, ZERO, psi(:,nb1), kdmx )
    ELSE
       CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, psi, &
                   kdmx, vc, nvecx, ZERO, psi(:,nb1), kdmx )
    ENDIF
     
    DO np = 1, notcnv
      psi(:,nbase+np) = -ew(nbase+np)*psi(:,nbase+np)
    ENDDO

    CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi, &
                kdmx, vc, nvecx, ONE, psi(1,nb1), kdmx )

    !write(*,*) 'psi(1,nb1) = ', psi(1,nb1)
    !stop 'ffr 252'

    ! approximate inverse iteration
    !CALL g_psi( npwx, npw, notcnv, npol, psi(1,1,nb1), ew(nb1) )
    
    ! "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
    ! order to improve numerical stability of subspace diagonalization
    ! (cdiaghg) ew is used as work array :
    !
    !      ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
    DO n = 1, notcnv
      nbn = nbase + n
      ew(n) = ddot( 2*npw, psi(:,nbn), 1, psi(:,nbn), 1 )
    ENDDO

    DO n = 1, notcnv
      psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
    ENDDO

    !write(*,*) 'psi(1,nb1) = ', psi(1,nb1)
    !stop 'ffr 272'

     
    ! here compute the hpsi and spsi of the new functions
    !CALL h_psi( npwx, npw, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
    !IF( uspp ) CALL s_psi( npwx, npw, notcnv, psi(1,1,nb1), spsi(1,1,nb1) )
    hpsi(:,nb1:nb1+notcnv-1) = matmul(H, psi(:,nb1:nb1+notcnv-1))
    if(uspp) spsi(:,nb1:nb1+notcnv-1) = matmul(S, psi(:,nb1:nb1+notcnv-1))

    ! update the reduced hamiltonian
    CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                kdmx, hpsi(:,nb1), kdmx, ZERO, hc(:,nb1), nvecx )
    


    IF( uspp ) THEN
      CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                  kdmx, spsi(:,nb1), kdmx, ZERO, sc(:,nb1), nvecx )
    ELSE
      CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                  kdmx, psi(:,nb1), kdmx, ZERO, sc(:,nb1), nvecx )
    ENDIF

    nbase = nbase + notcnv
    DO n = 1, nbase
      ! the diagonal of hc and sc must be strictly real 
      hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
      sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
      !
      DO m = n + 1, nbase
        hc(m,n) = CONJG( hc(n,m) )
        sc(m,n) = CONJG( sc(n,m) )
      ENDDO
    ENDDO

    write(100,*) real(Hc)
    write(101,*) imag(Hc)

    write(102,*) real(Sc)
    write(103,*) imag(Sc)

    write(*,*) 'Eigenvalues of reduced Hamiltonian before cdiaghg'
    do i=1,nbase
      write(*,'(1x,I4,F18.10)') i, ew(i)
    enddo

    ! diagonalize the reduced hamiltonian
    CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )

    write(*,*) 'nbase = ', nbase
    write(*,*) 'nvec = ', nvec
    write(*,*) 'Eigenvalues of reduced Hamiltonian after cdiaghg'
    do i=1,nbase
      write(*,'(1x,I4,F18.10)') i, ew(i)
    enddo

    write(104,*) real(vc)
    write(105,*) imag(vc)

    !return ! early return

    ! test for convergence
    WHERE( btype(1:nvec) == 1 )
      conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
    ELSEWHERE
      conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
    END WHERE

    write(*,*) 'ew[1:Nvec] = ', ew(1:nvec)
    write(*,*) 'e[1:Nvec] = ', e(1:nvec)

    notcnv = COUNT( .NOT. conv(:) )
    write(*,*) 'notcnv = ', notcnv
    write(*,*) 'conv = ', conv

    stop 'ffr 344'

    e(1:nvec) = ew(1:nvec)

    ! if overall convergence has been achieved, or the dimension of
    ! the reduced basis set is becoming too large, or in any case if
    ! we are at the last iteration refresh the basis set. i.e. replace
    ! the first nvec elements with the current estimate of the
    ! eigenvectors;  set the basis dimension to nvec.
    IF( notcnv == 0 .OR. &
        nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN

      CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
                   psi, kdmx, vc, nvecx, ZERO, evc, kdmx )

      IF( notcnv == 0 ) THEN
        ! all roots converged: return
        EXIT iterate
        !
      ELSEIF( dav_iter == maxter ) THEN
        ! ... last iteration, some roots not converged: return
        !
        !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
        !!!     &   " eigenvalues not converged")' ) notcnv
        !
        EXIT iterate
        !
      ENDIF

      ! refresh psi, H*psi and S*psi
      psi(:,1:nvec) = evc(:,1:nvec)
      IF( uspp ) THEN
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, spsi, &
                    kdmx, vc, nvecx, ZERO, psi(:,nvec+1), kdmx )
        spsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
      ENDIF

      CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, hpsi, &
                  kdmx, vc, nvecx, ZERO, psi(:,nvec+1), kdmx )

      hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)

      ! refresh the reduced hamiltonian 
      nbase = nvec
      hc(:,1:nbase) = ZERO
      sc(:,1:nbase) = ZERO
      vc(:,1:nbase) = ZERO
      DO n = 1, nbase
        hc(n,n) = CMPLX( e(n), 0.0_DP ,kind=DP)
        sc(n,n) = ONE
        vc(n,n) = ONE
      ENDDO

    ENDIF

  ENDDO iterate

  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  RETURN
  !
END SUBROUTINE 

