# !INPUT/OUTPUT PARAMETERS:
#   sol  : speed of light in atomic units (in,real)
#   kpa  : quantum number kappa (in,integer)
#   e    : energy (in,real)
#   nr   : number of radial mesh points (in,integer)
#   r    : radial mesh (in,real(nr))
#   vr   : potential on radial mesh (in,real(nr))
#   nn   : number of nodes (out,integer)
#   g0   : m th energy derivative of the major component multiplied by r
#          (out,real(nr))
#   g1   : radial derivative of g0 (out,real(nr))
#   f0   : m th energy derivative of the minor component multiplied by r
#          (out,real(nr))
#   f1   : radial derivative of f0 (out,real(nr))
function rdiracint!(sol, kpa, e, nr, r, vr, g0, g1, f0, f1)

    # rescaling limit
    rsc = 1.0e100

    @assert nr >= 4

    # inverse speed of light
    ci = 1.0/sol
    
    # electron rest energy
    e0 = sol^2
    t1 = 2.0*e0 + e
    
    # determine the r -> 0 boundary values of F and G
    t2 = kpa/r[1]
    t3 = ci*(t1 - vr[1])
    t4 = ci*(vr[1] - e)
    f0[1] = 1.0
    f1[1] = 0.0
    g0[1] = (f1[1] - t2*f0[1])/t4
    g1[1] = t3*f0[1] - t2*g0[1]

    # extrapolate to the first four points
    g1[2:4] .= g1[1]
    f1[2:4] .= f1[1]
    
    nn = 0 # number of nodes
    for ir in 2:nr
        t2 = kpa/r[ir]
        t3 = ci*(t1 - vr[ir])
        t4 = ci*(vr[ir] - e)
        ir0 = ir - 3
        if ir0 < 1
            ir0 = 1
        end
        @views g1[ir] = poly3(r[ir0:ir0+2], g1[ir0:ir0+2], r[ir])
        @views f1[ir] = poly3(r[ir0:ir0+2], f1[ir0:ir0+2], r[ir])
        # integrate to find wavefunction
        @views g0[ir] = poly4i( r[ir0:ir0+3], g1[ir0:ir0+3], r[ir]) + g0[ir0]
        @views f0[ir] = poly4i( r[ir0:ir0+3], f1[ir0:ir0+3], r[ir]) + f0[ir0]
        # compute the derivatives
        g1[ir] = t3*f0[ir] - t2*g0[ir]
        f1[ir] = t4*g0[ir] + t2*f0[ir]
        # integrate for correction
        @views g0[ir] = poly4i(r[ir0:ir0+3], g1[ir0:ir0+3], r[ir]) + g0[ir0]
        @views f0[ir] = poly4i(r[ir0:ir0+3], f1[ir0:ir0+3], r[ir]) + f0[ir0]
        # compute the derivatives again
        g1[ir] = t3*f0[ir] - t2*g0[ir]
        f1[ir] = t4*g0[ir] + t2*f0[ir]
        # check for overflow
        if ( (abs(g0[ir]) > rsc) || (abs(g1[ir]) > rsc) ||
             (abs(f0[ir]) > rsc) || (abs(f1[ir]) > rsc) )
            # set the remaining points and return
            g0[ir:nr] .= g0[ir]
            g1[ir:nr] .= g1[ir]
            f0[ir:nr] .= f0[ir]
            f1[ir:nr] .= f1[ir]
            return nn, e
        end
        # check for node
        if g0[ir-1]*g0[ir] < 0.0
            nn = nn + 1
        end
    end
    return nn, e
end
