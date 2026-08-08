includet("compute_phius.jl")

function debug_gener_pseudo_01(; NiterMax=100)

    ld1x_input = create_input_Si()
    #ld1x_input = create_input_Pd()

    Zval = ld1x_input.Zval
    Zed = ld1x_input.Zed
    Nspin = ld1x_input.Nspin
    Nwf = ld1x_input.Nwf
    nn = ld1x_input.nn
    ll = ld1x_input.ll
    oc = ld1x_input.oc
    core_state = ld1x_input.core_state

    Enl = zeros(Float64, Nwf)

    @assert ld1x_input.iswitch == 1

    rmax = 100.0
    xmin = -7.0
    dx = 1.25e-2
    ibound = false # default

    # Initialize radial grid
    grid = RadialGrid(rmax, Zval, xmin, dx, ibound)
    println("Nrmesh = ", grid.Nrmesh)

    Nrmesh = grid.Nrmesh
    v0 = zeros(Float64, Nrmesh)
    vxt = zeros(Float64, Nrmesh)
    Vpot = zeros(Float64, Nrmesh, Nspin)
    
    Rhoe = zeros(Float64, Nrmesh, Nspin)
    Rhoe_radial = zeros(Float64, Nrmesh, Nspin)
    V_h = zeros(Float64, Nrmesh)
    Vxc = zeros(Float64, Nrmesh, Nspin)
    VHxc = zeros(Float64, Nrmesh, Nspin)
    VHxc_new = zeros(Float64, Nrmesh, Nspin)
    epsxc = zeros(Float64, Nrmesh)
    Vnew = zeros(Float64, Nrmesh, Nspin)
    psi = zeros(Float64, Nrmesh, Nwf) # spin index dropped for the moment
    
    xc_calc = LibxcXCCalculator() # default using VWN
    ispin = 1
    
    starting_potential!(
        Nrmesh, Zval, Zed,
        Nwf, oc, nn, ll,
        grid.r, Enl, v0, vxt, Vpot
    )
    if Nspin == 2
        # XXX Same starting potential for spinpol case
        for i in 1:Nrmesh
           Vpot[i,2] = Vpot[i,1]
        end
    end
    # Define VHxc input
    for i in 1:Nrmesh
        VHxc[i,ispin] = Vpot[i,ispin] + Zed/grid.r[i]
    end


    # Solve for all states
    thresh0 = 1.0e-10
    nstop = 0
    mode = 1 # for lschps
    ze2 = -Zval # should be 2*Zval in Ry unit
    diff_V = Inf

    for iterSCF in 1:NiterMax

        println("\niterSCF = ", iterSCF)

        # FIXME: simplify this
        for iwf in 1:Nwf
            # Skip this iwf if occupation number is negative (unbound state)
            if oc[iwf] < 0.0
                Enl[iwf] = 0.0
                @views fill!(psi[:,iwf], 0.0)
                continue
                # The code below should not be executed 
            end
            @views psi1 = psi[:,iwf] # zeros wavefunction
            if ld1x_input.rel == 1
                Enl[iwf], nstop = lschps!( mode, Zval, thresh0, 
                    grid, nn[iwf], ll[iwf], Enl[iwf], Vpot, psi1
                )
            else   
                Enl[iwf], nstop = ascheq!(
                    nn[iwf], ll[iwf], Enl[iwf], grid, Vpot, ze2, thresh0, psi1, nstop
                )
            end
        end

        println("Energy levels:")
        for iwf in 1:Nwf
            @printf("%3d %18.10f\n", iwf, Enl[iwf])
        end

        #
        # calculate charge density (spherical approximation)
        #
        fill!(Rhoe, 0.0)
        for iwf in 1:Nwf, ir in 1:Nrmesh
            # this is for ispin=1
            Rhoe[ir,ispin] += oc[iwf] * psi[ir,iwf]^2
        end
        integRho = PWDFT.integ_simpson(Nrmesh, Rhoe[:,1], grid.rab) 
        println("integRho = ", integRho)

        radial_poisson_solve!(0, 2, grid, Rhoe, V_h)
        println("sum V_h = ", sum(V_h))

        Rhoe_radial[:] .= Rhoe[:] ./ grid.r2[:] ./ (4π) # 

        calc_epsxc_Vxc_VWN!(
            xc_calc, Rhoe_radial,
            epsxc,
            Vxc
        )
        println("sum Rhoe_radial = ", sum(Rhoe_radial))
        println("sum epsxc = ", sum(epsxc))
        println("sum Vxc = ", sum(Vxc))

        for i in 1:Nrmesh
            VHxc_new[i,ispin] = V_h[i] + Vxc[i,ispin]
            Vnew[i,ispin] = -Zed/grid.r[i] + vxt[i] + VHxc_new[i,ispin]
        end

        diff_V = LinearAlgebra.norm(VHxc .- VHxc_new)
        println("diff_V = ", diff_V)
        if diff_V < 1e-10
            println("!!!! CONVERGED !!!")
            break
        end

        for i in 1:Nrmesh
            # Mix
            VHxc[i,ispin] = 0.5*VHxc_new[i,ispin] + 0.5*VHxc[i,ispin]
            # Set new potential
            Vpot[i,ispin] = -Zed/grid.r[i] + VHxc[i,ispin]
        end

    end

    # We want to pseudize Vpot

    # rcloc is cutoff for 
    rcloc = ld1x_input.rcloc
    ir_loc = 0
    for i in 1:Nrmesh
        if grid.r[i] < rcloc
            ir_loc = i
        end
    end
    println("ir_loc = ", ir_loc)
    if ir_loc % 2 == 0
        println("Found ir_loc is even number, making it odd")
        ir_loc += 1
    end
    @assert ir_loc > 1
    @assert ir_loc < Nrmesh

    println("ir_loc = ", ir_loc)

    fae = Vpot[ir_loc]
    f1ae = deriv_7pts(Vpot, ir_loc, grid.r[ir_loc], grid.dx)
    f2ae = deriv2_7pts(Vpot, ir_loc, grid.r[ir_loc], grid.dx)

    log_der_ae = f1ae/fae
    ncn = 2
    𝓁 = 0
    flag = 0
    xc = zeros(Float64, 8)
    @views ld1x_find_qi!(grid, log_der_ae, xc[4:4+ncn-1], ir_loc, 𝓁, ncn, flag)
    
    j1 = zeros(Float64, Nrmesh, 8)
    norm_fact = zeros(Float64, 2)
    # compute the functions
    for ic in 1:2
        # CALL sph_bes(ik+1, grid%r, xc(3+nc), 0, j1(1,nc))
        for ir in 1:(ir_loc+1)
            j1[ir,ic] = sphericalbesselj(0, xc[3+ic]*grid.r[ir])
        end
        norm_fact[ic] = Vpot[ir_loc] / j1[ir_loc,ic]
        for ir in 1:(ir_loc+1)
            j1[ir,ic] = j1[ir,ic]*norm_fact[ic]
        end
    end
    #
    # compute the second derivative and impose continuity of zero, 
    # first and second derivative
    bm = zeros(Float64, 2) # the derivative of the bessel
    for ic in 1:2
        p1aep1 = ( j1[ir_loc+1,ic] - j1[ir_loc,ic] ) / ( grid.r[ir_loc+1] - grid.r[ir_loc] )
        p1aem1 = ( j1[ir_loc,ic] - j1[ir_loc-1,ic] ) / ( grid.r[ir_loc] - grid.r[ir_loc-1] )
        bm[ic] = (p1aep1 - p1aem1)*2 / ( grid.r[ir_loc+1] - grid.r[ir_loc-1] )
    end
  
    xc[2] = ( f2ae - bm[1] ) / ( bm[2] - bm[1] )
    xc[1] = 1.0 - xc[2]

    V_Ps_loc = zeros(Float64, Nrmesh)
    #
    # define the v_out function
    for ir in 1:ir_loc
        V_Ps_loc[ir] = xc[1]*j1[ir,1] + xc[2]*j1[ir,2]
    end
    
    # Beyond rcloc V_Ps_loc should be the same as Vpot
    for ir in (ir_loc+1):Nrmesh
        V_Ps_loc[ir] = Vpot[ir]
    end

    # We calculate rhoe core here
    #
    # calculates core charge density
    #
    rhov = zeros(Float64, Nrmesh)
    rhoc = zeros(Float64, Nrmesh)
    for ir in 1:grid.Nrmesh
        for iwf in 1:Nwf
            if ld1x_input.rel == 2
                # This is the case of full relativistic (Dirac equation)
                #XXX We need to reshape psi spinor
                if core_state[iwf]
                    rhoc[ir] += oc[iwf]*( psi[ir,1,iwf]^2 + psi[ir,2,iwf]^2 )
                else
                    rhov[ir] += oc[iwf]*( psi[ir,1,iwf]^2 + psi[ir,2,iwf]^2 )
                end
            else
                # scalar relativistic and non-relativistic
                if core_state[iwf]
                    rhoc[ir] += oc[iwf]*psi[ir,iwf]^2
                else
                    rhov[ir] += oc[iwf]*psi[ir,iwf]^2
                end
            end
        end
    end
    totrho = integ_0_inf_dr(rhoc, grid, Nrmesh, 2)
    println("totrho = ", totrho)

    # rcore determined by the condition  rhoc(rcore) = 2*rhov(rcore)
    ir_rcore = 0
    rcore = ld1x_input.rcore
    if ld1x_input.rcore > 0.0
        for ir in 1:Nrmesh
            if grid.r[ir] > rcore
                ir_rcore = ir
                break
            end
        end
    else 
        for ir in 1:Nrmesh
            if rhoc[ir] < 2.0*rhov[ir]
                println("at this ir=$ir rhoc=$(rhoc[ir]) rhov=$(rhov[ir]) 2*rhov=$(2*rhov[ir])")
                ir_rcore = ir
                break
            end
        end
    end
    #ir_rcore = 794
    println("ir_rcore = ", ir_rcore)
    rcore = grid.r[ir_rcore]
    println("rcore = ", rcore)

    # This is not used?
    dd1 = rhoc[ir_rcore+1]/grid.r2[ir_rcore+1] - rhoc[ir_rcore]/grid.r2[ir_rcore]
    drho = dd1/grid.dx/grid.r[ir_rcore]
    println("drho = ", drho)

    aeccharge = zeros(Float64, Nrmesh)
    aeccharge[:] = rhoc[:]
    #xc = zeros(Float64, 8) # again?
    compute_phius!(grid, 1, ir_rcore, aeccharge, rhoc, xc, 0, " ")

    totrho = integ_0_inf_dr(rhoc, grid, Nrmesh, 2)
    println("integ rhoc = ", totrho)

    @infiltrate

    return
end
