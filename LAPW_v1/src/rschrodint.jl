#=
! !INPUT/OUTPUT PARAMETERS:
!   sol : speed of light in atomic units (in,real)
!   l   : angular momentum quantum number (in,integer)
!   e   : energy (in,real)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   P₀  : m th energy derivative of P (out,real(nr))
!   P₁  : radial derivative of P₀ (out,real(nr))
!   Q₀  : m th energy derivative of Q (out,real(nr))
!   Q₁  : radial derivative of Q₀ (out,real(nr))
=#
function rschrodint!(
    l::Int64, E::Float64,
    r, Vr,
    P₀, P₁, Q₀, Q₁;
    sol=137.035999084
)
#=
  IMPLICIT NONE
  ! arguments
  REAL(8), INTENT(in) :: sol
  INTEGER, INTENT(in) :: l
  REAL(8), INTENT(in) :: e
  INTEGER, INTENT(in) :: nr
  REAL(8), INTENT(in) :: r(nr),vr(nr)
  INTEGER, INTENT(out) :: nn
  REAL(8), INTENT(out) :: P₀(nr),P₁(nr)
  REAL(8), INTENT(out) :: Q₀(nr),Q₁(nr)
  ! local variables
  INTEGER :: ir,ir0
  REAL(8) :: ri,t1,t2,t3,t4
=#
    t1 = 1.0/sol^2
    t2 = l*(l+1)
    # determine the r -> 0 boundary values of P and Q
    ri = 1.0/r[1]
    t3 = 2.0 + t1*( E - Vr[1] )
    t4 = t2/(t3*r[1]^2) + Vr[1] - E
    Q₀[1] = 1.0
    Q₁[1] = 0.0
    P₀[1] = ( Q₁[1] + Q₀[1]*ri)/t4
    P₁[1] = t3*Q₀[1] + P₀[1]*ri
    # extrapolate to the first four points
    P₁[2:4] .= P₁[1]
    Q₁[2:4] .= Q₁[1]
    #
    nn = 0
    nr = size(r, 1)
    for ir in 2:nr
        ri = 1.0/r[ir]
        t3 = 2.0 + t1*( E - Vr[ir] )
        t4 = t2/(t3*r[ir]^2) + Vr[ir] - E
        ir0 = ir - 3
        if ir0 < 1
            ir0 = 1
        end
        @views P₁[ir] = poly3( r[ir0:ir0+2], P₁[ir0:ir0+2], r[ir] )
        @views Q₁[ir] = poly3( r[ir0:ir0+2], Q₁[ir0:ir0+2], r[ir] )
        # integrate to find wavefunction
        @views P₀[ir] = poly4i( r[ir0:ir0+3], P₁[ir0:ir0+3], r[ir] ) + P₀[ir0]
        @views Q₀[ir] = poly4i( r[ir0:ir0+3], Q₁[ir0:ir0+3], r[ir] ) + Q₀[ir0]
        # compute the derivatives
        P₁[ir] = t3*Q₀[ir] + P₀[ir]*ri
        Q₁[ir] = t4*P₀[ir] - Q₀[ir]*ri
        # integrate for correction
        @views P₀[ir] = poly4i( r[ir0:ir0+3], P₁[ir0:ir0+3], r[ir] ) + P₀[ir0]
        @views Q₀[ir] = poly4i( r[ir0:ir0+3], Q₁[ir0:ir0+3], r[ir] ) + Q₀[ir0]
        # compute the derivatives again
        P₁[ir] = t3*Q₀[ir] + P₀[ir]*ri
        Q₁[ir] = t4*P₀[ir] - Q₀[ir]*ri
        # check for overflow
        if ( abs(P₀[ir]) > 1.0e100 ) || ( abs(P₁[ir]) > 1.0e100) ||
           ( abs(Q₀[ir]) > 1.0e100 ) || ( abs(Q₁[ir]) > 1.0e100) 
            P₀[ir:nr] = P₀[ir]
            P₁[ir:nr] = P₁[ir]
            Q₀[ir:nr] = Q₀[ir]
            Q₁[ir:nr] = Q₁[ir]
            return nn, E
        end
        # check for node
        if P₀[ir-1]*P₀[ir] < 0.0
            nn += 1
        end
    end

    return nn, E
end
