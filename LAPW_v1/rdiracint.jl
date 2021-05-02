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

function poly3(xa, ya, x)
    # arguments
    # real(8) xa(3),ya(3),x
    # local variables

    # evaluate the polynomial coefficients
    x0 = xa[1]
    x1 = xa[2] - x0
    x2 = xa[3] - x0
    y0 = ya[1]
    y1 = ya[2] - y0
    y2 = ya[3] - y0
    t0 = 1.0/(x1*x2*(x2-x1))
    t1 = x1*y2
    t2 = x2*y1
    c1 = x2*t2 - x1*t1
    c2 = t1 - t2
    t1 = x - x0
    
    # evaluate the polynomial
    return y0 + t0*t1*(c1 + c2*t1)
end

function poly4i(xa, ya, x)
    # arguments
    #real(8), intent(in) :: xa(4),ya(4),x

    # evaluate the polynomial coefficients
    x0 = xa[1]
    x1 = xa[2] - x0
    x2 = xa[3] - x0
    x3 = xa[4] - x0
    y0 = ya[1]
    y1 = ya[2] - y0
    y2 = ya[3] - y0
    y3 = ya[4] - y0
    t4 = x1 - x2
    t5 = x1 - x3
    t6 = x2 - x3
    t1 = x1*x2*y3
    t2 = x2*x3*y1
    t3 = x1*x3
    t0 = 1.0/(x2*t3*t4*t5*t6)
    t3 = t3*y2
    c3 = t1*t4 + t2*t6 - t3*t5
    t4 = x1^2
    t5 = x2^2
    t6 = x3^2
    c2 = t1*(t5 - t4) + t2*(t6 - t5) + t3*(t4 - t6)
    c1 = t1*(x2*t4 - x1*t5) + t2*(x3*t5 - x2*t6) + t3*(x1*t6 - x3*t4)
    t1 = x - x0
    # integrate the polynomial
    return t1*(y0 + t0*t1*(0.5*c1 + t1*(c2/3 + 0.25*c3*t1) ) )
end

