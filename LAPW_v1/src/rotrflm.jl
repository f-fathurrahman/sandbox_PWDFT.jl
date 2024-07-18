function rotrflm!(R, lmax, n, ld, rflm1, rflm2)
    #=
    ! !INPUT/OUTPUT PARAMETERS:
!   rot   : rotation matrix (in,real(3,3))
!   lmax  : maximum angular momentum (in,integer)
!   n     : number of functions to rotate (in,integer)
!   ld    : leading dimension (in,integer)
!   rflm1 : coefficients of the real spherical harmonic expansion for each
!           function (in,real(ld,n))
!   rflm2 : coefficients of rotated functions (out,complex(ld,n))
    =#

    #=
    ! arguments
  REAL(8), intent(in) :: rot(3,3)
  INTEGER, intent(in) :: lmax,n,ld
  REAL(8), intent(in) :: rflm1(ld,*)
  REAL(8), intent(out) :: rflm2(ld,*)
  ! local variables
  INTEGER l,lm,nm,p
  REAL(8) det,rotp(3,3)
  REAL(8) ang(3),angi(3)
  ! automatic arrays
  REAL(8) d(ld,ld)
    =#

    @assert lmax >= 0
    @assert n > 0
    if n == 0
        return
    end

    ld = size(rflm1, 1)

    ang = zeros(Float64, 3, 3)
    angi = zeros(Float64, 3, 3)
    d = zeros(Float64, ld, ld)

    detR = det(R)

    # make the rotation proper
    p = 1
    if detR < 0.0
        p = -1
    end
    pR = p*R
    
    # compute the Euler angles of the rotation matrix
    roteuler!(pR, ang)
    # inverse rotation: the function is to be rotated, not the spherical harmonics
    angi[1] = -ang[3]
    angi[2] = -ang[2]
    angi[3] = -ang[1]
    # determine the rotation matrix for real spherical harmonics
    #CALL rlmrot(p, angi, lmax, ld, d)
    rlmrot!(p, angi, lmax, d)
    # apply rotation matrix
    for l in 0:lmax
        nm = 2*l + 1
        lm = l^2 + 1
        #CALL dgemm('N', 'N', nm, n, nm, 1.0, d(lm,lm), ld, rflm1(lm,1), ld, 0.d0, rflm2(lm,1), ld)
        rflm2[lm:(lm+nm-1), 1:nm] = d[lm:(lm+nm-1), lm:(lm+n-1)] * rflm1[lm:(lm+n-1), 1:nm]
        # (nm x n) (n x nm) -> (nm,nm) 
        #
    end
    return
end