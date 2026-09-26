function ascheqps_Ry!(
    nam, ℓ, E0, grid, Vpot, y, beta_prj, ddd, qq, lls, idx_rbeta;
    TOL = 1e-11, NmaxIter = 20
)
    # TOL might be too small ?
#=
numerical integration of a generalized radial schroedinger equation,
using Numerov with outward and inward integration and matching.
Works for both norm-conserving nonlocal and US pseudopotentials
Requires in input a good estimate "E0" of the energy
=#

    println("\nEnter ascheqps")
    println("n=$nam, ℓ=$ℓ, input E0 = $E0")

    Nrmesh = grid.Nrmesh
    Nbeta = size(beta_prj, 2)
    @assert Nbeta == size(ddd, 1)
    @assert Nbeta == size(ddd, 2)

    # set up constants and allocate variables the 
    fun = zeros(Float64, Nrmesh)
    f = zeros(Float64, Nrmesh)
    el = zeros(Float64, Nrmesh)
    c = zeros(Float64, Nrmesh)
    work = zeros(Float64, Nbeta)
    
    nstop = 0
    ir_start = 0
    E = E0
    eup = 0.3*E
    elw = 1.3*E
    ndcr = nam - ℓ - 1

    ddx12 = grid.dx^2/12.0
    nst = 2*(ℓ + 1)
    sqlhf = (ℓ + 0.5)^2
    #
    # series developement of the potential near the origin
    for ir in 1:4
        y[ir] = Vpot[ir]
    end
    #
    println()
    println("Before radial_grid_series")
    @printf("y[1] = %18.10f\n", y[1])
    @printf("y[2] = %18.10f\n", y[2])
    @printf("y[3] = %18.10f\n", y[3])
    @printf("y[4] = %18.10f\n", y[4])
    println()
    #
    b = zeros(Float64, 4) # originally b(0:3)
    radial_grid_series!( y, grid.r, grid.r2, b )
    #
    println()
    @printf("b[1] = %18.10f\n", b[1])
    @printf("b[2] = %18.10f\n", b[2])
    @printf("b[3] = %18.10f\n", b[3])
    @printf("b[4] = %18.10f\n", b[4])
    println()
    #
    #  set up the f-function and determine the position of its last
    #  change of sign
    #  f < 0 (approximatively) means classically allowed   region
    #  f > 0         "           "        "      forbidden   "
    #
    for iterSch in 1:NmaxIter
        println("starting iterSch=$iterSch, elw=$elw, E=$E, eup=$eup")
        idx_r = 1
        #f[1] = ddx12*( 2 * grid.r2[1] * (Vpot[1] - E) + sqlhf ) # XXX change to Ha
        f[1] = ddx12*( grid.r2[1] * (Vpot[1] - E) + sqlhf ) # XXX This is Ry
        for ir in 2:Nrmesh
            #f[ir] = ddx12*( 2 * grid.r2[ir] * (Vpot[ir] - E) + sqlhf ) # XXX change to Ha
            f[ir] = ddx12*( grid.r2[ir] * (Vpot[ir] - E) + sqlhf ) # XXX This is in Ry
            if ( f[ir] != abs(f[ir])*sign(f[ir-1]) ) && (ir < Nrmesh-5)
                idx_r = ir
            end
        end
        if (idx_r == 1) || (grid.r[idx_r] > 4.0)
            idx_r = round(Int64, Nrmesh*3/4)
        end
        #
        if idx_r >= Nrmesh-2
            error("No point found for matching")
            # probably something wrong with the potential
        end
        #
        # determine if idx_r is sufficiently large
        for ibeta in 1:Nbeta
            if (lls[ibeta] == ℓ) && (idx_rbeta[ibeta] > idx_r)
                idx_r = idx_rbeta[ibeta] + 3
            end
        end
        println("idx_r = ", idx_r)
        #
        # if everything is ok continue the integration and define f
        for ir in 1:Nrmesh
            f[ir] = 1.0 - f[ir]
        end
        #
        # determination of the wave-function in the first two points by
        # series developement
        #
        # no coulomb divergence in the origin for a pseudopotential
        ze2 = 0.0 
        #start_scheq_Ha!( ℓ, E, b, grid, ze2, y )
        start_scheq_Ry!( ℓ, E, b, grid, ze2, y )
        #
        # outward integration before idx_r
        integrate_outward_Ry!( ℓ, E, grid, f, b, y, beta_prj, ddd, qq, lls, idx_rbeta, idx_r)
        ncross = 0
        ymx = 0.0
        for ir in 2:(idx_r-1)
            # XXX why not simply compare the sign?
            if y[ir] != abs(y[ir])*sign(y[ir+1]) # XXX this is float comparison
                ncross += 1
            end
            ymx = max(ymx, abs(y[ir+1]))
        end
        # If at this point the number of nodes is wrong it means that something
        # is probably wrong in the calling routines. A ghost might be present
        # in the pseudopotential. With a nonlocal pseudopotential there is no
        # node theorem so strictly speaking the following instructions are
        # wrong but sometimes they help so we keep them here.
        if ndcr != ncross # `first` should be inout argument
            println("Warning: n=$nam l=$ℓ, expecting $ndcr nodes, found $ncross")
            println("Setting wfc to zero for this iteration")
        end
        #
        if ndcr < ncross
            # too many crossings. E is an upper bound to the true eigenvalue.
            # increase abs(E)
            #
            eup = E
            E = 0.9*elw + 0.1*eup
            println("too many crossing: ncross=$ncross, ndcr=$ndcr")
            fill!(y, 0.0)
            ymx = 0.0
            @goto LABEL300
        elseif ndcr > ncross
            #
            # too few crossings. E is a lower bound to the true eigenvalue.
            # decrease abs(E)
            #
            elw = E
            E = 0.9*eup + 0.1*elw
            println("Too few crossing ncross=$ncross, ndcr=$ndcr")
            fill!(y, 0.0) #XXX Need this? 
            ymx = 0.0
            @goto LABEL300
        end
        #
        # inward integration up to idx_r
        #
        ir_start = integrate_inward!(grid, f, y, c, el, idx_r)
        #
        # if necessary, improve the trial eigenvalue by the cooley's procedure.
        # jw cooley math of comp 15,363(1961)
        #
        fe = (12.0 - 10.0*f[idx_r])*y[idx_r] - f[idx_r-1]*y[idx_r-1] - f[idx_r+1]*y[idx_r+1]
        #
        # adjust the normalization if needed
        if ymx >= 1.0e10
            @views y[:] = y[:]/ymx
        end
        #
        # calculate the normalization
        for ibeta in 1:Nbeta
            if (ℓ == lls[ibeta]) # also need to check jj for relativistic case
                idx_r_l = idx_rbeta[ibeta]
                for ir in 1:idx_r_l
                    fun[ir] = beta_prj[ir,ibeta]*y[ir]*sqrt(grid.r[ir])
                end
                work[ibeta] = integ_0_inf_dr(fun, grid, idx_r_l, nst)
            else
                work[ibeta] = 0.0
            end
        end
        #
        for ir in 1:ir_start
            fun[ir] = y[ir]*y[ir]*grid.r[ir]
        end
        ss = integ_0_inf_dr(fun, grid, ir_start, nst)
        for ibeta in 1:Nbeta, jbeta in 1:Nbeta
            ss += qq[ibeta,jbeta]*work[ibeta]*work[jbeta]
        end
        dfe = -y[idx_r]*f[idx_r]/grid.dx/ss
        de = -fe*dfe #  in Ry
        #de = -fe*dfe/2 # Hartree?
        epsE = abs(de/E)
        println("iterSch = $iterSch E=$E de=$de")
        if abs(de) < TOL
            println("CONVERGED: at iterSch = $iterSch E = $E de = $de")
            @goto LABEL600
        end
        #
        if epsE > 0.25
            de = 0.25*de/epsE
        end
        #
        if de > 0.0
            elw = E
        end
        if de < 0.0
            eup = E
        end
        E = E + de
        #
        if E > eup
            E = 0.9*eup + 0.1*elw
        end
        #
        if E < elw
            E = 0.9*elw + 0.1*eup
        end
        @label LABEL300
    end
    nstop = 1
    
    if ir_start == 0
        @goto LABEL900 # return?
    end
  
    @label LABEL600
    #  
    # exponential tail of the solution if it was not computed
    #
    if ir_start < Nrmesh
        for ir in ir_start:Nrmesh-1
            if y[ir] == 0.0 # XXX Float comparison
                y[ir+1] = 0.0
            else
                yln = log(abs(y[ir]))
                xp = -sqrt(12.0*abs(1.0-f[ir]))
                expn = yln + xp
                if expn < -80.0
                    y[ir+1] = 0.0
                else
                    y[ir+1] = abs(exp(expn))*sign(y[ir])
                end
            end
        end
    end
    #
    # normalize the eigenfunction as if they were norm conserving. 
    # If this is a US PP the correct normalization is done outside this routine.
    #
    for ir in 1:Nrmesh
        el[ir] = grid.r[ir]*y[ir]*y[ir]
    end
    ss = integ_0_inf_dr(el, grid, Nrmesh, nst)
    if ss > 0.0
        ss = sqrt(ss)
        for ir in 1:Nrmesh
            y[ir] = sqrt(grid.r[ir]) * y[ir] / ss
        end
        E0 = E
    else
        nstop = 1
    end
    @label LABEL900
    return E

end


