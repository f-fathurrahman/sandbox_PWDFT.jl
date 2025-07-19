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
        norm_u = dot(u.^2, dr) # use radial_wavefunc_norm
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

# Finite-difference derivative of radial function.
function radial_derivative!(rgrid, n_g, dndr_g)
    #
    npts = rgrid.npts
    #
    dndr_g[1] = n_g[2] - n_g[1] # forward difference
    #
    # centered difference
    for i in 2:npts-1
        dndr_g[i] = 0.5*(n_g[i+1] - n_g[i-1])
    end
    dndr_g[npts] = n_g[npts] - n_g[npts-1] # backward difference
    @. dndr_g /= rgrid.dr
    return
end

function radial_wavefunc_norm(u, dr)
    idx = abs.(u) .>= 1e-160
    return dot( u[idx].^2, dr[idx] ) 
end


# FIXME: some arguments can be obtained from rgrid
function solve_radial_sch!(
    rgrid, N, r, dr, vr, d2gdr2, n_j, l_j,e_j, u_j;
    scalarrel=true
)
    c2 = @. -(r/dr)^2 # beware this can become a matrix without @.
    c10 = @. -d2gdr2 * r^2  # first part of c1 vector

    if scalarrel
        r2dvdr = zeros(Float64, N)
        radial_derivative!(rgrid, vr, r2dvdr) # calc dvdr
        r2dvdr[:] .*= r
        r2dvdr[:] .-= vr
    else
        r2dvdr = nothing
    end
    #println("sum abs r2dvdr = ", sum(abs.(r2dvdr)))
    # XXX: r2dvdr is local to this function

    # solve for each quantum state separately
    Nstates = length(n_j)
    #for ist in 1:Nstates
    for ist in 4:4
        n = n_j[ist]
        l = l_j[ist]
        E = e_j[ist]
        u = u_j[ist]
        #
        nodes = n - l - 1  # analytically expected number of nodes
        Δ = -0.2*E
        nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
        #return
        # adjust eigenenergy until u has the correct number of nodes
        while nn != nodes
            dif_nn = sign(nn - nodes)
            while dif_nn == sign(nn - nodes)
                print("Different sign: integrate again")
                E -= dif_nn * Δ
                nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
            end
            Δ = Δ/2
        end
        println("No of nodes already good")
        println("nodes should be ", nodes)

        # normalize
        norm_u = radial_wavefunc_norm(u, dr) # no check of abs u
        u[:] .*= 1.0 / sqrt(norm_u)

        #return # early return for debug

        println("\n---- Adjusting eigenenergy until smooth ----")

        # adjust eigenenergy until u is smooth at the turning point
        dE = 1.0
        iterNo = 1
        MAX_ITER_ADJ = 50
        while abs(dE) > 1e-9
            println("\niterNo=$iterNo E = $E dE = $dE")
            norm_u = radial_wavefunc_norm(u, dr)
            #
            println("norm_u = ", norm_u)
            u[:] .*= 1.0 / sqrt(norm_u)
            println("norm_u after = ", radial_wavefunc_norm(u, dr))
            dE = 0.5 * A / norm_u
            x = abs(dE/E)
            println("x = ", x)
            if x > 0.1
                dE = dE * 0.1 / x # 0.1 * E/dE * dE
                println("dE is adjusted to = ", dE)
            end
            E -= dE # update
            @assert E < 0.0
            nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
            if nn != nodes
                error("nn=$nn should be equal to $nodes")
            end
            iterNo += 1
            if iterNo >= MAX_ITER_ADJ
                println("!!!!!! Too much iteration")
                break
            end
        end
        e_j[ist] = E # Set
        norm_u = radial_wavefunc_norm(u, dr)
        u[:] .*= 1.0 / sqrt(norm_u) # normalize
    end
    return
end


function integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel; gmax=nothing)

    println("\nENTER integ_radial_sch")

    # FIXME: r2dvdr can be nothing
    if scalarrel
        α = 1 / 137.036  # fine-structure constant
        x = 0.5 * α^2  # x = 1 / (2c^2)
        Mr = r * (1.0 + x * E) - x * vr
    else
        Mr = r
    end
    #
    c0 = @. l * (l + 1) + 2 * Mr * (vr - E * r)
    if isnothing(gmax) && all(c0 .> 0)
        error("Bad initial electron density guess!")
    end

    c1 = c10
    if scalarrel
        @. c0 += x * r2dvdr / Mr
        @. c1 = c10 - x * r * r2dvdr / (Mr * dr)
    end
    #@printf("sum c10 = %18.10e\n", sum(c10))
    #@printf("sum c0 = %18.10e\n", sum(c0))
    #@printf("sum c1 = %18.10e\n", sum(c1))
    #@printf("sum c2 = %18.10e\n", sum(c2))


    # vectors needed for numeric integration of diff. equation
    fm = @. 0.5*c1 - c2
    fp = @. 0.5*c1 + c2
    f0 = @. c0 - 2*c2

    #@printf("sum fm = %18.10e\n", sum(fm))
    #@printf("sum fp = %18.10e\n", sum(fp))
    #@printf("sum f0 = %18.10e\n", sum(f0))    

    if isnothing(gmax)
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[end] = 1.0
        # perform backwards integration from infinity to the turning point
        g = length(u) - 1
        u[end-1] = f0[end]*u[end] / fm[end]
        while c0[g] > 0.0  # this defines the classical turning point
            u[g-1] = ( f0[g]*u[g] + fp[g]*u[g+1] ) / fm[g]
            if u[g-1] < 0.0
                println("WARNING: There should't be a node here!  Use a more negative")
                return 100, nothing
            end
            if u[g-1] > 1e100
                @info "Scaling u in integ_radial_sch"
                u[:] .*= 1e-100
            end
            g -= 1
        end
        # stored values of the wavefunction and the first derivative
        # at the turning point
        gtp = g + 1 # ???
        utp = u[gtp] # ????
        if gtp == length(u) # some errors or warnings ?
            @info "Early return in integ_radial_sch 237"
            return 100, 0
        end
        dudrplus = 0.5*( u[gtp+1] - u[gtp-1] ) / dr[gtp]
    else
        gtp = gmax
    end
    println("integ_radial_sch: After backward integration gtp=$gtp, c0=$(c0[gtp]) u=$(u[gtp])")

    # set boundary conditions at r -> 0
    u[1] = 0.0
    u[2] = 1.0

    # perform forward integration from zero to the turning point
    g = 2
    nodes = 0
    # integrate one step further than gtp
    # (such that dudr is defined in gtp)
    while g <= gtp+1
        u[g+1] = (fm[g]*u[g-1] - f0[g]*u[g]) / fp[g]
        if u[g+1] * u[g] < 0
            nodes += 1
        end
        g += 1
    end
    # what is use case of this?
    if !isnothing(gmax)
        return #???
    end
    println("my_eigen_shoot: after inward integration u=$(u[gtp])")

    # scale first part of wavefunction, such that it is continuous at gtp
    u[1:gtp+2] .*= utp / u[gtp] # utp from backward integ

    println("my_eigen_shoot: after scale u=$(u[gtp])")

    # determine size of the derivative discontinuity at gtp
    dudrminus = 0.5 * (u[gtp+1] - u[gtp-1]) / dr[gtp]
    A = (dudrplus - dudrminus) * utp
    
    println("nodes = ", nodes)
    @printf("utp = %18.10e\n", utp)
    @printf("dudrplus = %18.10e\n", dudrplus)
    @printf("dudrminus = %18.10e\n", dudrminus)
    @printf("A = %18.10e\n", A)

    println("EXIT integ_radial_sch\n")

    return nodes, A
end