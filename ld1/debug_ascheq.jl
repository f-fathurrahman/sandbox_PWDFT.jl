using Printf

include("RadialGrid.jl")
include("start_scheq.jl")

import PyPlot
const plt = PyPlot

mutable struct AtomicElectronsConfig
    Zval::Float64
    Nspin::Int64
    Nwf::Int64
    # These quantities might depend on Nspin,
    # but we assume that Nspin=1 for the moment
    nn::Vector{Int64}
    ll::Vector{Int64}
    Focc::Vector{Int64}
    Enl::Vector{Float64}
end

function init_atom_Si()
    Zval = 14.0
    Nspin = 1
    Nwf = 5
    nn = [1, 2, 2, 3, 3]
    ll = [0, 0, 1, 0, 1] 
    Focc = [2.0, 2.0, 6.0, 2.0, 2.0]
    Enl = zeros(Float64, Nwf)

    # Make sure that we don't mess up with nn, ll, and Focc (to some extent)
    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(Focc) == Nwf

    return AtomicElectronsConfig(
        Zval, Nspin, Nwf, nn, ll, Focc, Enl
    )

end

function init_radial_grid( Zval::Float64 )
    rmax = 100.0
    xmin = -7.0 # iswitch = 1
    dx = 0.008 # iswitch = 1
    ibound = false # default
    # Initialize radial grid
    grid = RadialGrid(rmax, Zval, xmin, dx, ibound)
    return grid
end


function starting_potential!(
  elec_config::AtomicElectronsConfig,
  grid::RadialGrid,
  V0, Vxt, Vpot;
  frozen_core=false, noscf=false
)
    
    # V0, Vxt are not used (?)
    # enne is not used (?)

    enne = 0.0
    # zz = max(Zed, Zval)
    # zed equal to zval (?_
    Zval = elec_config.Zval
    zz = Zval
    nn = elec_config.nn
    ll = elec_config.ll
    Focc = elec_config.Focc
    Enl = elec_config.Enl
    Nwf = elec_config.Nwf
    Nspin = elec_config.Nspin

    Nrmesh = grid.Nrmesh
    r = grid.r

    for iwf in 1:Nwf
       oce = max(0.0, Focc[iwf])
       enne = enne + oce
       zen = 0.0
       for jwf in 1:Nwf
            oce = max(0.0, Focc[jwf])
            if nn[jwf] < nn[iwf]
                zen = zen + oce
            end
            if (nn[jwf] == nn[iwf]) && (ll[jwf] <= ll[iwf])
                zen = zen + oce
            end
        end
        zen = max(zz - zen + 1.0, 1.0)
        if (abs(Enl[iwf]) < 1.e-7) || (!frozen_core)
            Enl[iwf] = -0.5*( zen/nn[iwf] )^2  # Ha unit
        end
    end
    
    for i in 1:Nrmesh
        # might have external potential, we set it to zero for the moment
        Vxt[i] = 0.0
        #
        x = r[i]*enne^(1.0/3.0)/0.885
        t = zz/(1.0 + sqrt(x)*(0.02747 - x*(0.1486 - 0.007298*x)) + x*(1.243 + x*(0.2302 + 0.006944*x)))
        t = max(1.0, t)
        V0[i] = -Zval/r[i] # Zed => Zval
        if noscf
            Vpot[i,1] = V0[i] + Vxt[i]
        else
            Vpot[i,1] = -t/r[i] + Vxt[i] # V0 is not used
        end
    end
    
    if Nspin == 2
        # Just use the same potential as up spin
        for i in 1:Nrmesh
          Vpot[i,2] = Vpot[i,1]
        end
    end
    return enne # not used?
end


function debug_ascheq!(
    Zval, 𝓃::Int64, 𝓁::Int64, Enl_guess::Float64,
    grid::RadialGrid,
    Vpot, y, thresh0;
    ispin=1
)


    Nrmesh = grid.Nrmesh
    println("input Enl_guess = ", Enl_guess)

    #
    # set up constants and initialize
    #
    c = zeros(Float64, Nrmesh)
    f = zeros(Float64, Nrmesh)
    el = zeros(Float64, Nrmesh)

    thresh = thresh0
    if Enl_guess < -5e2
        thresh = thresh0*10.0
    end
    
    ddx12 = grid.dx^2 / 12.0
    l1 = 𝓁 + 1
    sqlhf = 0.5*(𝓁 + 0.5)^2 # Ha unit
    ndcr = 𝓃 - 𝓁 - 1
    
    Enl = Enl_guess
    # set initial lower and upper bounds to the eigenvalue
    Enl_up = Vpot[Nrmesh] + sqlhf/grid.r2[Nrmesh]
    Enl_lw = Enl_up
    for i in 1:Nrmesh
       Enl_lw = min( Enl_lw, Vpot[i,ispin] + sqlhf/grid.r2[i] )
    end
    
    nstop = 200
    if Enl_up == Enl_lw
        error("Should go to 900")
    end
    
    if Enl_lw > Enl_up
        Enl = 0.9*Enl_up + 0.1*Enl_lw
    end
    
    if Enl < Enl_lw
        Enl = 0.9*Enl_lw + 0.1*Enl_up
    end

    println()
    println("Initial lower and upper bounds to the eigenvalue (in Ha)")
    println()
    @printf("Enl_guess = %18.10f\n", Enl)
    @printf("Enl_up    = %18.10f\n", Enl_up)
    @printf("Enl_lw    = %18.10f\n", Enl_lw)
    
    #
    #  series developement of the potential near the origin
    #
    for i in 1:4
        y[i] = Vpot[i] - Zval/grid.r[i]
    end
    b = zeros(Float64,4) # XXX: b is originally b(0:3)
    radial_grid_series!( y, grid.r, grid.r2, b )

    println()
    println("Near origin: r, Vpot, yi, b (in Ha)")
    println()
    for i in 1:4
        @printf("%18.10e %18.10e %18.10e %18.10e\n", grid.r[i], Vpot[i], y[i], b[i])
    end


    NmaxIter = 50

    # needed outside the loop
    nstart = 0
    f2 = 0.0

    println()

    # iterSch = 1
    for iterSch in 1:NmaxIter
    
        println("===============================")
        println("iterSch = ", iterSch)
        println("===============================")

        nstop = 300
    
        # set up the f-function and determine the position of its last
        # change of sign
        #
        # f < 0 (approximatively) means classically allowed   region
        # f > 0         "           "        "      forbidden   "
        #
        #
        # Using Numerov algorithm
        #
        ik = 0
        f[1] = ddx12*( grid.r2[1]*( Vpot[1] - Enl ) + sqlhf)
        for i in 2:Nrmesh
            f[i] = ddx12*( grid.r2[i]*( Vpot[i] - Enl ) + sqlhf)
            if f[i] != abs(f[i])*sign(f[i-1])
                ik = i
            end
        end
        println("f = ", f[1:4])
        println("ik = ", ik)
        nstop = 302
    
        # XXX: What's this?
        if ik >= (Nrmesh - 2)
            error("ik is too big: " * string(ik))
        end

        for i in 1:Nrmesh
            f[i] = 1.0 - 2*f[i] # convert to Ry?
        end
        for i in 1:4
            @printf("f = %18.10f\n", f[i])
        end
        fill!(y, 0.0)
  
        # Determination of the wave-function in the first two points
        # by series developement
        # TODO: add reference(s)
        xl1 = 𝓁 + 1.0
        x4l6 = 4.0*𝓁 + 6.0
        b0e = b[1] - Enl # Ha unit
        c1 = Zval/xl1  # in Ha?
        c2 = (c1*Zval + b0e)/x4l6 # Ha

        start_scheq!( 𝓁, Enl, b, grid, Zval, y )

        @printf("After start_scheq!\n")
        @printf("y[1] = %18.10f\n", y[1])
        @printf("y[2] = %18.10f\n", y[2])

        # start outward integration and count number of crossings
        ncross = 0
        ymx = 0.0
        println("\nStarting outward integration")
        println("ik = ", ik)
        println("fn = ", f[1:4])
        for i in 2:(ik-1)
            # Numerov algorithm
            y[i+1] = ( ( 12.0 - 10.0*f[i] )*y[i] - f[i-1]*y[i-1] ) / f[i+1]
            #
            # Check for crossing here
            #
            if y[i] != abs(y[i])*sign(y[i+1])
                println("Found crossing")
                ncross = ncross + 1
            end
            # Track maximum value of y (wavefunction)
            ymx = max( ymx, abs(y[i+1]) )
        end
        
        plt.clf()
        plt.plot(grid.r[1:ik-1], y[1:ik-1], color="blue", label="outward")
        
        @printf("ymx = %18.10f\n", ymx)
        println("ncross = ", ncross)

        # matching radius has been reached going out. if ncross is not
        # equal to ndcr, modify the trial eigenvalue.
        if ndcr < ncross
            # Too many crossings.
            # Enl is an upper bound to the true eigenvalue.
            # Increase abs(Enl)
            println("Too many crossing")
            Enl_up = Enl
            rap = ( Float64(ncross + l1)/𝓃 )^2
            Enl = (Enl - Vpot[Nrmesh] )*rap + Vpot[Nrmesh]
            if Enl < Enl_lw
                Enl = 0.9*Enl_lw + 0.1*Enl_up
            end
            # continue # need this?
        elseif ndcr > ncross
            # too few crossings.
            # Enl is a lower bound to the true eigenvalue.
            # Decrease abs(e)
            println("Too few crossing")
            elw = Enl
            rap = ( Float64(ncross + l1)/𝓃 )^2
            Enl = ( Enl - Vpot[Nrmesh] )*rap + Vpot[Nrmesh]
            if Enl > eup
                Enl = 0.9*Enl_up + 0.1*Enl_lw
            end
            # go to 300
            # continue
        end
    
        println("ndcr = ", ndcr)

        # Prepare inward integration
        # charlotte froese can j phys 41,1895(1963)
        #
        # Start at min(rmax, 10*rmatch)
        #
        nstart = Nrmesh
        ns = 10
        rstart = ns*grid.r[ik]
        if rstart < grid.r[Nrmesh]
            for i in ik:Nrmesh
                nstart = i
                if grid.r[i] >= rstart
                    #go to 403
                    break
                end
            end
            # 403 
            nstart = round(Int64, floor(nstart/2))
            nstart = 2*nstart + 1
        end
        println("nstart = ", nstart)
        #if iter == 1
        #    println("exit ffr 231")
        #    exit()
        #end


        # set up a, l, and c vectors
        n = ik + 1
        el[n] = 10.0*f[n] - 12.0
        c[n] = -f[ik]*y[ik]
        n2 = ik + 2
        for n in n2:nstart
           di = 10.0*f[n] - 12.0
           el[n] = di - f[n]*f[n-1]/el[n-1]
           c[n] = -c[n-1]*f[n-1]/el[n-1]
        end
        println("el[nstart] = ", el[nstart])
        println("c[nstart] = ", c[nstart])
        #if iter == 1
        #    println("exit ffr 247")
        #    exit()
        #end

        # Start inward integration by the Froese's tail procedure
        expn = exp( -sqrt( 12.0*abs(1.0 - f[nstart-1]) ) )
        y[nstart-1] = c[nstart-1]/( el[nstart-1] + f[nstart]*expn )
        y[nstart] = expn*y[nstart-1]
        for n in range(nstart-2, stop=ik+1, step=-1) #nstart-2,ik+1,-1
            y[n] = ( c[n] - f[n+1]*y[n+1])/el[n]
        end
        @printf("y = %18.10f\n", y[ik+1])

        plt.plot(grid.r[ik+1:nstart], y[ik+1:nstart], color="red", label="inward")
        plt.legend()
        plt.grid(true)
        plt.title("Enl = " * string(Enl))
        plt.savefig("IMG_trial_psi_" * string(iterSch) * ".png", dpi=150)


        # If necessary, improve the trial eigenvalue by the cooley's procedure.
        # J.W. Cooley, Math of Comp 15, 363 (1961)
        fe = ( 12.0 - 10.0*f[ik] )*y[ik] - f[ik-1]*y[ik-1] - f[ik+1]*y[ik+1]
        @printf("fe = %18.10f\n", fe)

        #  calculate the normalization
        if ymx >= 1.0e10
           for i in 1:Nrmesh
              y[i] = y[i]/ymx
           end
        end
        @printf("ymx = %18.10f\n", ymx)        
        @printf("y[1] = %18.10f\n", y[1])
        @printf("y[Nrmesh] = %18.10f\n", y[Nrmesh])

        plt.clf()
        plt.plot(grid.r[1:nstart], y[1:nstart], label="trial psi (normalized)")
        plt.legend()
        plt.grid(true)
        plt.savefig("IMG_trial_psi_normalized_" * string(iterSch) * ".png", dpi=150)
    
        a0 = 1.0/(2*𝓁 + 3)
        a1 = c1/(𝓁 + 2)
        a2 = (c1*c1 + c2 + c2)/(2*𝓁 + 5)
        sum0 = (a0 + grid.r[1]*( a1 + grid.r[1]*a2 ) )*grid.r[1]^(2*𝓁 + 3)
        nst2 = nstart - 2
        f2 = grid.r2[1]*y[1]*y[1]
        ss = grid.r[1]*f2/(2*l1 + 1)
        # Integrate using Simpson's rule?
        for n in range(1, stop=nst2, step=2)
            f0 = f2
            f1 = grid.r2[n+1]*y[n+1]*y[n+1]
            f2 = grid.r2[n+2]*y[n+2]*y[n+2]
            ss = ss + f0 + f2 + 4.0*f1
        end
        @printf("ss = %18.10f\n", ss)
        ss = sum0 + grid.dx*ss/3.0
        dfe = -y[ik]*f[ik]/grid.dx/ss
        de = -fe*dfe*0.5 # Ha unit?
        eeps = abs(de/Enl)
        
        @printf("ss = %18.10f\n", ss)
        @printf("de = %18.10f\n", de)

        println("iterSch = ", iterSch, " Enl = ", Enl,  " de = ", de)
        if abs(de) < 2*thresh
            println("*** CONVERGED ***")
            break
        end
    
        if eeps > 0.25
            println("Updating de in 341")
            de = 0.25*de/eeps
        end

        if de > 0.0
            println("Updating Enl_lw")
            Enl_lw = Enl
        end
    
        if de < 0.0
            println("Updating Enl_up")
            Enl_up = Enl
        end
  
        Enl = Enl + de
        println("de  = ", de)
        println("Enl = ", Enl)
  
        if Enl > Enl_up
            Enl = 0.9*Enl_up + 0.1*Enl_lw
        end
    
        if Enl < Enl_lw
            Enl = 0.9*Enl_lw + 0.1*Enl_up
        end
  
        nstop = 50 # ???
    
        println("New data:")
        @printf("Enl    = %18.10f\n", Enl)
        @printf("Enl_up = %18.10f\n", Enl_up)
        @printf("Enl_lw = %18.10f\n", Enl_lw)

    end


    println("Final Enl = ", Enl)

    # Set the values for tail
    for i in nstart:(Nrmesh-1)
        y[i+1] = 0.0
        if y[i] == 0.0 # FIXME
            continue
        end        
        yln = log(abs(y[i]))
        xp = -sqrt(12.0*abs(1.0 - f[i]))
        expn = yln + xp
        if expn < -80.0
            continue
        end
        y[i+1] = abs(exp(expn))*sign(y[i])
    end

    # normalize the eigenfunction and exit
    sum1 = 0.0
    println("nstart for normalization: ", nstart)
    for i in range(nstart, stop=Nrmesh-2, step=2)
        f0 = f2
        f1 = grid.r2[i+1] * y[i+1] * y[i+1]
        f2 = grid.r2[i+2] * y[i+2] * y[i+2]
        sum1 = sum1 + f0 + f2 + 4.0*f1
    end
    ss = ss + grid.dx*sum1/3.0
    ss = sqrt(ss)
    for i in 1:Nrmesh
        y[i] = grid.sqrtr[i]*y[i]/ss
    end

    # Check integration
    xl1 = 𝓁 + 1.0
    x4l6 = 4.0*𝓁 + 6.0
    b0e = b[1] - Enl # Ha unit
    c1 = Zval/xl1
    c2 = (c1*Zval + b0e)/x4l6 # Ha
    a0 = 1.0/(2*𝓁 + 3)
    a1 = c1/(𝓁 + 2)
    a2 = (c1*c1 + c2 + c2)/(2*𝓁 + 5)
    sum0 = (a0 + grid.r[1]*( a1 + grid.r[1]*a2 ) )*grid.r[1]^(2*𝓁 + 3)
    nst2 = nstart - 2
    f2 = grid.r2[1]*y[1]*y[1]
    ss = grid.r[1]*f2/(2*l1 + 1)
    # Integrate using Simpson's rule?
    for i in range(1, stop=Nrmesh-2, step=2)
        f0 = f2
        f1 = grid.r2[i+1]*y[i+1]*y[i+1]
        f2 = grid.r2[i+2]*y[i+2]*y[i+2]
        ss = ss + f0 + f2 + 4.0*f1
    end
    @printf("ss = %18.10f\n", ss)
    ss = sum0 + grid.dx*ss/3.0
    println("Check norm: ", ss) # ???

    
    if nstop < 100
        println("Should go to 900")
        return Enl, nstop
    end
    
    nstop = 0
    return Enl, nstop
end



function test_main()
    elec_config = init_atom_Si()
    grid = init_radial_grid(elec_config.Zval)
    
    Nrmesh = grid.Nrmesh
    V0 = zeros(Float64, Nrmesh)
    Vxt = zeros(Float64, Nrmesh)
    Vpot = zeros(Float64, Nrmesh, 2)
    enne = starting_potential!(
        elec_config, grid,
        V0, Vxt, Vpot
    )
    println("enne = ", enne) # should be equal to Zval
    println("elec_config.Enl = ")
    println(elec_config.Enl)

    #plt.clf()
    #plt.plot(grid.r, Vpot[:,1], label="Vpot")
    #plt.xlim(0.0, 0.2)
    #plt.savefig("IMG_Vpot.png", dpi=150)

    # Solve for all states
    ze2 = -elec_config.Zval # using negative sign and should be 2*Zval in Ry unit
    thresh0 = 1.0e-10
    Nwf = elec_config.Nwf
    psi = zeros(Float64, Nrmesh, Nwf)
    iwf = 2
    #for iwf in 1:Nwf
        @views psi1 = psi[:,iwf] # zeros wavefunction
        elec_config.Enl[iwf], nstop = debug_ascheq!(
            ze2,
            elec_config.nn[iwf],
            elec_config.ll[iwf],
            elec_config.Enl[iwf],
            grid, Vpot,
            psi1, thresh0
        )
    #end

end

test_main()


