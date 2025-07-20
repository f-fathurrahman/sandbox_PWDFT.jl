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
    for ist in 1:Nstates

        #println("\nSolving radial Schroedinger equation for ist = $ist")

        n = n_j[ist]
        l = l_j[ist]
        E = e_j[ist]
        u = u_j[ist]
        #
        # initial radial integration
        nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
        #
        #println("--- Early return in solve_radial_sch for debug !!!!")
        #return
        #
        # adjust eigenenergy until u has the correct number of nodes
        nodes = n - l - 1  # analytically expected number of nodes
        Δ = -0.2*E # XXX Is this safe?
        #println("nn=$nn nodes=$nodes")
        iterNode = 0
        NiterNodeMax = 100
        while nn != nodes
            dif_nn = sign(nn - nodes)
            #println("dif_nn = ", dif_nn)
            # as long as no difference in diff_nn
            # dif_nn will be positive if for current E we will have more nodes
            while dif_nn == sign(nn - nodes)
                E -= dif_nn * Δ
                #println("Same sign but nn != nodes")
                #println("updating E to $(E)")
                nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
                #println("nn=$nn nodes=$nodes")
            end
            #println("Different sign or already match: nn=$nn nodes=$nodes")
            Δ *= 1/2 # reduce
            #
            iterNode += 1
            if iterNode >= NiterNodeMax
                println("WARNING: Could not converge no of nodes")
                break
            end
        end
        #println("No of nodes already good")
        #println("nodes should be $nodes nn=$nn")

        # normalize
        norm_u = radial_wavefunc_norm(u, dr) # no check of abs u
        u[:] .*= 1.0 / sqrt(norm_u)
        #println("norm_u after normalization = ", norm_u)

        #return # early return for debug

        #println("\n---- Adjusting eigenenergy until smooth ----")

        # adjust eigenenergy until u is smooth at the turning point
        dE = 1.0
        iterNo = 1
        MAX_ITER_ADJ = 100
        is_converged = false
        for iterNo in 1:MAX_ITER_ADJ
            #println("\niterNo=$iterNo Current E = $E dE = $dE")
            norm_u = radial_wavefunc_norm(u, dr)
            #
            #println("norm_u = ", norm_u)
            u[:] .*= 1.0 / sqrt(norm_u)
            #println("norm_u after normalization = ", radial_wavefunc_norm(u, dr))
            #
            dE = 0.5 * A/norm_u
            x = abs(dE/E)
            #println("dE=$dE x=$x")
            # why this?
            if x > 0.1
                dE = dE * 0.1 / x # 0.1 * E/dE * dE
                #println("dE is adjusted to = $dE because c > 0.1")
            end
            E -= dE # update
            #println("E is updated to $E")
            @assert E < 0.0
            nn, A = integ_radial_sch!(u, l, vr, E, r2dvdr, r, dr, c10, c2, scalarrel)
            # nn should not changed
            if nn != nodes
                error("nn=$nn should be equal to $nodes")
            end
            iterNo += 1
            if abs(dE) <= 1e-9
                #println("CONVERGED!!!")
                is_converged = true
                break 
            end
        end
        #println("is_converged = ", is_converged)
        if !is_converged
            println("WARNING: Cannot converge adjutst E")
        end
        e_j[ist] = E # Set
        norm_u = radial_wavefunc_norm(u, dr)
        u[:] .*= 1.0 / sqrt(norm_u) # normalize
    end
    return
end