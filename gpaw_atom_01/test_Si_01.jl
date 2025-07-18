using Infiltrator

using LinearAlgebra 
using Printf

using PWDFT: LibxcXCCalculator, calc_epsxc_Vxc_LDA!

includet("AERadialGrid.jl")

function radial_integrate(rgrid, a_xg; n=0)
    @assert n >= -2
    r = rgrid.r
    dr = rgrid.dr
    # exclude the first radial point?
    return dot( a_xg[2:end], (r.^(2 + n) .* dr)[2:end] ) * (4π)
end

# Initialize with Slater function
function init_wavefunc!(atsymb, r, dr, l_j, E_j, u_j)
    #
    Nstates = length(l_j)
    special_case = ["Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au"]
    NradialPoints = length(r)
    #
    for ist in 1:Nstates
        l = l_j[ist]
        E = E_j[ist]
        u = u_j[ist]
        if atsymb in special_case
            a = sqrt(-4.0*E)
        else
            a = sqrt(-2.0*E)
        end
        for i in 1:NradialPoints
            u[i] = r[i]^(1 + l) * exp(-a * r[i])
        end
        norm_u = dot(u.^2, dr)
        u[:] .*= 1.0 / sqrt(norm_u)
    end
    return
end

# Return the electron charge density divided by 4 pi
function calc_density!(r, f_j, u_j, rhoe)
    NradialPoints = length(rhoe)
    Nstates = length(f_j)
    pref = 1/(4π)
    for ist in 1:Nstates, i in 2:NradialPoints
        if abs(u_j[ist][i]) >= 1e-160
            rhoe[i] += pref * f_j[ist] * u_j[ist][i]^2 / r[i]^2
        end
    end
    rhoe[1] = rhoe[2]
    return
end

# input n*r*dr, output vr=vh*r
function solve_radial_hartree!(
    l::Int64, nrdr, r, vr
)

    NradialPoints = length(r)    
    p = 0.0
    q = 0.0
    fill!(vr, 0.0)
    for i in range(NradialPoints, stop=1, step=-1)
        rl = r[i]^l
        dp = nrdr[i]/rl
        rlp1 = rl * r[i]
        dq = nrdr[i] * rlp1
        vr[i] = (p + 0.5*dp)*rlp1 - (q + 0.5*dq)/rl
        p += dp
        q += dq
    end

    #println("p = ", p, " q = ", q)
    
    vr[1] = 0.0
    f = 4π / (2*l + 1)
    #for (int g = 1; g < M; g++)
    for i in 2:NradialPoints
        vr[i] = f * (vr[i] + q / r[i]^l)
    end

    return
end

function solve_radial_sch(
    rgrid, N, r, dr, vr, d2gdr2, n_j, l_j,e_j, u_j;
    scalarrel=true
)
    c2 = -(r/dr)^2
    c10 = -d2gdr2 * r^2  # first part of c1 vector

    if scalarrel
        r2dvdr = np.zeros(N)
        rgd.derivative(vr, r2dvdr)
        r2dvdr *= r
        r2dvdr -= vr
    else
        r2dvdr = nothing
    end
    # XXX: r2dvdr is local to this function

    # solve for each quantum state separately
    Nstates = length(n_j)
    for ist in 1:Nstates
        n = n_j[ist]
        l = l_j[ist]
        E = e_j[ist]
        u = u_j[ist]
        #
        nodes = n - l - 1  # analytically expected number of nodes
        Δ = -0.2 * E
        nn, A = integ_radial_sch(u, l, vr, e, r2dvdr, r, dr, c10, c2, scalarrel)
        # adjust eigenenergy until u has the correct number of nodes
        while nn != nodes
            dif_nn = sign(nn - nodes)
            while dif_nn == sign(nn - nodes)
                e -= dif_nn * Δ
                nn, A = integ_radial_sch(u, l, vr, e, r2dvdr, r, dr, c10, c2, scalarrel)
            end
            Δ = Δ/2
        end

        # adjust eigenenergy until u is smooth at the turning point
        dE = 1.0
        while abs(dE) > 1e-9
            norm_u = dot(u.^2, dr) # no check of abs u
            u .*= 1.0 / sqrt(norm_u)
            de = 0.5 * A / norm_u
            x = abs(dE/E)
            if x > 0.1
                dE *= 0.1 / x
            end
            E = E - dE
            @assert E < 0.0
            nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
        end
        e_j[ist] = E # Set
        u .*= 1.0 / sqrt( dot(u.^2, dr) ) # normalize
    end
    return
end


function integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel; gmax=nothing)
    if scalarrel
        x = 0.5 * α^2  # x = 1 / (2c^2)
        Mr = r * (1.0 + x * E) - x * vr
    else
        Mr = r
    end
    c0 = l * (l + 1) + 2 * Mr * (vr - e * r)
    if isnothing(gmax) && all(c0 .> 0)
        error("Bad initial electron density guess!")
    end

    c1 = c10
    if scalarrel
        c0 += x * r2dvdr / Mr
        c1 = c10 - x * r * r2dvdr / (Mr * dr)
    end

    # vectors needed for numeric integration of diff. equation
    fm = 0.5 * c1 - c2
    fp = 0.5 * c1 + c2
    f0 = c0 - 2 * c2

    if isnothing(gmax)
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[-1] = 1.0

        # perform backwards integration from infinity to the turning point
        g = len(u) - 2
        u[-2] = u[-1] * f0[-1] / fm[-1]
        while c0[g] > 0.0  # this defines the classical turning point
            u[g - 1] = (f0[g] * u[g] + fp[g] * u[g + 1]) / fm[g]
            if u[g - 1] < 0.0
                prinlnt("WARNING: There should't be a node here!  Use a more negative")
                return 100, nothing
            end
            if u[g - 1] > 1e100
                u *= 1e-100
            end
            g -= 1
        end

        # stored values of the wavefunction and the first derivative
        # at the turning point
        gtp = g + 1
        utp = u[gtp]
        if gtp == len(u) - 1
            return 100, 0.0
        end
        dudrplus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    else
        gtp = gmax
    end

    # set boundary conditions at r -> 0
    u[1] = 0.0
    u[2] = 1.0

    # perform forward integration from zero to the turning point
    g = 1
    nodes = 0
    # integrate one step further than gtp
    # (such that dudr is defined in gtp)
    while g <= gtp
        u[g + 1] = (fm[g] * u[g - 1] - f0[g] * u[g]) / fp[g]
        if u[g + 1] * u[g] < 0
            nodes += 1
        end
        g += 1
    end
    if !isnothing(gmax)
        return
    end

    # scale first part of wavefunction, such that it is continuous at gtp
    u[:gtp + 2] *= utp / u[gtp]

    # determine size of the derivative discontinuity at gtp
    dudrminus = 0.5 * (u[gtp + 1] - u[gtp - 1]) / dr[gtp]
    A = (dudrplus - dudrminus) * utp

    return nodes, A
end




function main_test_Si_01()

    atsymb = "Si"
    Z = 14

    # Input configuration
    n_j =  [1, 2, 2, 3, 3] # main/principal quantum number
    l_j =  [0, 0, 1, 0, 1] # orbital
    f_j =  [2, 2, 6, 2, 2] # occupation
    e_j =  [-65.184426, -5.075056, -3.514938, -0.398139, -0.153293]

    maxnodes = maximum([n - l - 1 for (n, l) in zip(n_j, l_j)])
    #print("maxnodes = ", maxnodes)

    gpernode = 150 # default?
    β = 0.4 # parameter for radial grid
    NradialPoints = gpernode*(maxnodes + 1)

    a = β/NradialPoints
    b = 1/NradialPoints
    ae_grid = AERadialGrid(a, b, npts=NradialPoints)

    r = ae_grid.r
    dr = ae_grid.dr
    d2gdr2 = ae_grid.d2gdr2

    Nstates = length(n_j)
    # Radial wave functions multiplied by radius:
    u_j = Vector{Vector{Float64}}(undef, Nstates)
    for ist in 1:Nstates
        u_j[ist] = zeros(Float64, NradialPoints)
    end

    # Effective potential multiplied by radius:
    vr = zeros(Float64, NradialPoints)

    # Electron density:
    rhoe = zeros(Float64, NradialPoints)

    vHr = zeros(Float64, NradialPoints)
    vXC = zeros(Float64, NradialPoints)
    epsxc = zeros(Float64, NradialPoints)

    # Initialize starting wavefunctions and calculate density from them
    init_wavefunc!(atsymb, r, dr, l_j, e_j, u_j)
    calc_density!(r, f_j, u_j, rhoe)

    solve_radial_hartree!(0, rhoe .* r .* dr, r, vHr)

    # add potential from nuclear point charge (v = -Z / r)
    vHr .-= Z # vHr is vH times r, so the potential to added is (-Z/r)*r => -Z

    xc_calc = LibxcXCCalculator(x_id=1, c_id=12)

    fill!(vXC, 0.0)
    calc_epsxc_Vxc_LDA!(xc_calc, rhoe, epsxc, vXC)

    Exc = radial_integrate(ae_grid, epsxc .* rhoe)
    @. vr = (vHr + vXC*r)

    @exfiltrate
end
