function rotrflm!(R, lmax, n, ld, rflm1_, rflm2_)
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

    println("ld = ", ld)

    rflm1 = reshape(rflm1_, ld, :)
    rflm2 = reshape(rflm2_, ld, :)

    ang = zeros(Float64, 3)
    angi = zeros(Float64, 3)
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
        idx1 = lm:(lm+nm-1)
        idx2 = 1:n
        idx3 = lm:(lm+nm-1)
        @info "nm = $nm lm = $lm"
        @info "idx1 = $(idx1)"
        @info "idx2 = $(idx2)"
        @info "idx3 = $(idx3)"
        #=
        C <- A * B
        On entry,  M  specifies  the number  of rows  of the  matrix
        op( A )  and of the  matrix  C.  M  must  be at least  zero.
        On entry,  N  specifies the number  of columns of the matrix
        op( B ) and the number of columns of the matrix C. N must be
        at least zero.
        On entry,  K  specifies  the number of columns of the matrix
        op( A ) and the number of rows of the matrix op( B ). K must
        be at least  zero.
        =#
        #CALL dgemm('N', 'N', nm, n, nm, 1.0, d(lm,lm), ld, rflm1(lm,1), ld, 0.d0, rflm2(lm,1), ld)
        rflm2[idx1, idx2] = d[idx1, idx3] * rflm1[idx3, idx2]
        # (nm x n) (n x nm) -> (nm,nm) 
        #
    end
    return
end