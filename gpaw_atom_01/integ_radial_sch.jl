function integ_radial_sch!(
    u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel::Bool;
    gmax=nothing
)
    # FIXME: r2dvdr can be nothing

    println("\nENTER integ_radial_sch with E=$E l=$l")

    # Calculate some arrays that depend on E
    if scalarrel
        α = 1 / 137.036  # fine-structure constant
        x = 0.5*α^2  # x = 1 / (2c^2)
        Mr = r*(1.0 + x*E) - x*vr
    else
        Mr = r
    end
    #
    #NradialPoints = length(r)
    #c0 = zeros(Float64, NradialPoints)
    #for i in 1:NradialPoints
    #    c0[i] = l*(l + 1) + 2*Mr[i]*(vr[i] - E*r[i])
    #end
    c0 = @. l * (l + 1) + 2 * Mr * (vr - E * r)
    if isnothing(gmax) && all(c0 .> 0)
        error("Bad initial electron density guess!")
    end

    c1 = copy(c10)
    if scalarrel
        @. c0 += x * r2dvdr / Mr
        @. c1 = c10 - x * r * r2dvdr / (Mr * dr)
    end
    @printf("sum c10 = %18.10e\n", sum(c10))
    @printf("sum c0 = %18.10e\n", sum(c0))
    @printf("sum c1 = %18.10e\n", sum(c1))
    @printf("sum c2 = %18.10e\n", sum(c2))


    # vectors needed for numeric integration of diff. equation
    fm = @. 0.5*c1 - c2
    fp = @. 0.5*c1 + c2
    f0 = @. c0 - 2*c2

    @printf("sum fm = %18.10e\n", sum(fm))
    @printf("sum fp = %18.10e\n", sum(fp))
    @printf("sum f0 = %18.10e\n", sum(f0))    

    if isnothing(gmax)
        # set boundary conditions at r -> oo (u(oo) = 0 is implicit)
        u[end] = 1.0
        # perform backwards integration from infinity to the turning point
        g = length(u) - 1
        u[end-1] = f0[end]*u[end] / fm[end] # because u[end+1] is zero
        while c0[g] > 0.0  # this defines the classical turning point
            u[g-1] = ( f0[g]*u[g] + fp[g]*u[g+1] ) / fm[g]
            if u[g-1] < 0.0
                println("WARNING: There should't be a node here!  Use a more negative")
                error()
                return 100, nothing
            end
            if u[g-1] > 1e100
                @info "Scaling u in integ_radial_sch"
                error()
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
    nodes = 0
    # integrate one step further than gtp
    # (such that dudr is defined in gtp)
    for g in 2:gtp+1
        u[g+1] = (fm[g]*u[g-1] - f0[g]*u[g]) / fp[g]
        # Counts number of nodes
        if u[g+1]*u[g] < 0.0 && (g <= gtp) 
            println("Found nodes between: $(g+1) and $(g)")
            println("Values: $(u[g+1]) and $(u[g])")
            nodes += 1
        end
    end
    # what is use case of this?
    if !isnothing(gmax)
        return #???
    end
    println("my_eigen_shoot: after inward integration u=$(u[gtp])")

    # scale first part of wavefunction, such that it is continuous at gtp
    u[1:gtp+2] .*= (utp/u[gtp]) # utp from backward integ

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