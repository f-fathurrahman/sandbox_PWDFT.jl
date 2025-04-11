

function test_scf_01(; NiterMax=100)

    #ld1x_input = create_input_Si()
    ld1x_input = create_input_Pd()

    Zval = ld1x_input.Zval
    Zed = ld1x_input.Zed
    Nspin = ld1x_input.Nspin
    Nwf = ld1x_input.Nwf
    nn = ld1x_input.nn
    ll = ld1x_input.ll
    oc = ld1x_input.oc

    Enl = zeros(Float64, Nwf)

    @assert ld1x_input.iswitch == 1

    rmax = 100.0
    if ld1x_input.rel == 1
        xmin = -8.0 # iswitch=1, rel=1
    else
        xmin = -7.0 # iswitch=1, rel=0
    end
    dx = 0.008 # iswitch=1, rel=1
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

    #@infiltrate

    return

end
