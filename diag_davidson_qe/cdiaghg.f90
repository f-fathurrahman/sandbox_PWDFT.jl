!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both ZHEGV and ZHEGVX

  !
  IMPLICIT NONE
  !
  integer, parameter :: dp=8
  
  ! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: n

  ! number of eigenstates to be calculate
  integer, intent(in) :: m

  ! leading dimension of h, as declared in the calling pgm unit  
  integer, intent(in) :: ldh


  ! actually intent(in) but compilers don't know and complain
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n)  ! matrix to be diagonalized
  complex(dp), intent(inout) :: s(ldh,n)  ! overlap matrix

  REAL(DP), INTENT(OUT) :: e(n) ! eigenvalues

  COMPLEX(DP), INTENT(OUT) :: v(ldh,m) ! eigenvectors (column-wise)
  
  !
  ! Local variables
  !
  ! mm = number of calculated eigenvectors
  INTEGER :: lwork, nb, mm, info, i, j

  REAL(DP) :: abstol
  
  ! various work space
  INTEGER, ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
  
  LOGICAL :: all_eigenvalues

! ILAENV returns optimal block size "nb"
  INTEGER, EXTERNAL :: ILAENV

  !
  ! ... save the diagonal of input S (it will be overwritten)
  !
  ALLOCATE( sdiag( n ) )
  DO i = 1, n
     sdiag(i) = DBLE( s(i,i) )
  END DO

  all_eigenvalues = ( m == n )
  write(*,*) 'all_eigenvalues = ', all_eigenvalues

  ! check for optimal block size
  nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
  IF ( nb < 1 .OR. nb >= n) THEN
    lwork = 2*n
  ELSE
    lwork = ( nb + 1 )*n
  ENDIF
  
  ALLOCATE( work( lwork ) )
  
  IF( all_eigenvalues ) THEN
    ALLOCATE( rwork( 3*n - 2 ) )
    ! calculate all eigenvalues (overwritten to v)
    v(:,:) = h(:,:)
    CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
  ELSE
    ALLOCATE( rwork( 7*n ) )
    ! save the diagonal of input H (it will be overwritten)
    ALLOCATE( hdiag( n ) )
    DO i = 1, n
      hdiag(i) = DBLE( h(i,i) )
    ENDDO
    ALLOCATE( iwork(5*n) )
    ALLOCATE( ifail(n) )

    ! ... calculate only m lowest eigenvalues
    abstol = 0.D0

    CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                 0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                 work, lwork, rwork, iwork, ifail, info )

    DEALLOCATE( ifail )
    DEALLOCATE( iwork )

    ! restore input H matrix from saved diagonal and lower triangle
    DO i = 1, n
      h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
      DO j = i + 1, n
        h(i,j) = CONJG( h(j,i) )
      END DO
      DO j = n + 1, ldh
        h(j,i) = ( 0.0_DP, 0.0_DP )
      END DO
    END DO
    DEALLOCATE( hdiag )
  
  ENDIF
  
  DEALLOCATE( rwork )
  DEALLOCATE( work )

  IF( info > n ) THEN
    write(*,*) 'info = ', info
    stop 'cdiaghg: S matrix not positive definite' ! ABS( info ) )
  ELSEIF( info > 0 ) THEN
    write(*,*) 'info = ', info
    stop 'cdiaghg: eigenvectors failed to converge' ! ABS( info ) )
  ELSEIF( info < 0 ) THEN
    stop 'cdiaghg: incorrect call to ZHEGV*' ! ABS( info ) )
  ENDIF

  ! restore input S matrix from saved diagonal and lower triangle
  DO i = 1, n
    s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
    DO j = i + 1, n
      s(i,j) = CONJG( s(j,i) )
    ENDDO
    DO j = n + 1, ldh
      s(j,i) = ( 0.0_DP, 0.0_DP )
    ENDDO
  ENDDO
  !
  DEALLOCATE( sdiag )
  RETURN

END SUBROUTINE cdiaghg


