# !INPUT/OUTPUT PARAMETERS:
#   sol  : speed of light in atomic units (in,real)
#   κ  : quantum number kappa (in,integer)
#   e    : energy (in,real)
#   nr   : number of radial mesh points (in,integer)
#   r    : radial mesh (in,real(nr))
#   vr   : potential on radial mesh (in,real(nr))
#   nn   : number of nodes (out,integer)
#   G₀   : m th energy derivative of the major component multiplied by r
#          (out,real(nr))
#   G₁   : radial derivative of G₀ (out,real(nr))
#   F₀   : m th energy derivative of the minor component multiplied by r
#          (out,real(nr))
#   F₁   : radial derivative of F₀ (out,real(nr))
function rdiracint!(
    κ::Int64, E::Float64,
    r, Vr,
    G₀, G₁, F₀, F₁;
    sol=137.035999084
)

    nr = size(r, 1)

    # rescaling limit
    rsc = 1.0e100

    @assert nr >= 4

    # inverse speed of light
    ci = 1.0/sol
    
    # electron rest energy
    e0 = sol^2
    t1 = 2*e0 + E
    
    # determine the r -> 0 boundary values of F and G
    t2 = κ/r[1]
    t3 = ci*( t1 - Vr[1] )
    t4 = ci*( Vr[1] - E )
    F₀[1] = 1.0
    F₁[1] = 0.0
    G₀[1] = ( F₁[1] - t2*F₀[1] )/t4
    G₁[1] = t3*F₀[1] - t2*G₀[1]

    # extrapolate to the first four points
    @views G₁[2:4] .= G₁[1]
    @views F₁[2:4] .= F₁[1]
    
    nn = 0 # number of nodes
    for ir in 2:nr
        t2 = κ/r[ir]
        t3 = ci*(t1 - Vr[ir])
        t4 = ci*(Vr[ir] - E)
        ir0 = ir - 3
        if ir0 < 1
            ir0 = 1
        end
        @views G₁[ir] = poly3(r[ir0:ir0+2], G₁[ir0:ir0+2], r[ir])
        @views F₁[ir] = poly3(r[ir0:ir0+2], F₁[ir0:ir0+2], r[ir])
        # integrate to find wavefunction
        @views G₀[ir] = poly4i(r[ir0:ir0+3], G₁[ir0:ir0+3], r[ir]) + G₀[ir0]
        @views F₀[ir] = poly4i(r[ir0:ir0+3], F₁[ir0:ir0+3], r[ir]) + F₀[ir0]
        # compute the derivatives
        G₁[ir] = t3*F₀[ir] - t2*G₀[ir]
        F₁[ir] = t4*G₀[ir] + t2*F₀[ir]
        # integrate for correction
        @views G₀[ir] = poly4i(r[ir0:ir0+3], G₁[ir0:ir0+3], r[ir]) + G₀[ir0]
        @views F₀[ir] = poly4i(r[ir0:ir0+3], F₁[ir0:ir0+3], r[ir]) + F₀[ir0]
        # compute the derivatives again
        G₁[ir] = t3*F₀[ir] - t2*G₀[ir]
        F₁[ir] = t4*G₀[ir] + t2*F₀[ir]
        # check for overflow
        if (abs(G₀[ir]) > rsc) || (abs(G₁[ir]) > rsc) ||
           (abs(F₀[ir]) > rsc) || (abs(F₁[ir]) > rsc)
            # set the remaining points and return
            @views G₀[ir:nr] .= G₀[ir]
            @views G₁[ir:nr] .= G₁[ir]
            #
            @views F₀[ir:nr] .= F₀[ir]
            @views F₁[ir:nr] .= F₁[ir]
            return nn
        end
        # check for node
        if G₀[ir-1]*G₀[ir] < 0.0
            nn += 1
        end
    end
    return nn
end
