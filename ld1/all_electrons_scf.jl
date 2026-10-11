function all_electrons_scf!(
    grid, Zval, Zed, Nspin,
    Nwf, oc, nn, ll, Enl,
    xc_calc,
    psi, Rhoe, Rhoe_radial, Vpot, rel;
    NiterMax = 100
)
    
    Nrmesh = grid.Nrmesh

    v0 = zeros(Float64, Nrmesh)
    vxt = zeros(Float64, Nrmesh)
    V_h = zeros(Float64, Nrmesh)
    V_xc = zeros(Float64, Nrmesh, Nspin)
    V_hxc = zeros(Float64, Nrmesh, Nspin)
    V_hxc_new = zeros(Float64, Nrmesh, Nspin)
    epsxc = zeros(Float64, Nrmesh)
    V_new = zeros(Float64, Nrmesh, Nspin)

    starting_potential!(
        Nrmesh, Zval, Zed,
        Nwf, oc, nn, ll,
        grid.r, Enl, v0, vxt, Vpot
    )

    ispin = 1 # XXX HARCODED
    # XXX: This is not
    if Nspin == 2
        # XXX Same starting potential for spinpol case
        for i in 1:Nrmesh
           Vpot[i,2] = Vpot[i,1]
        end
    end

    # Define V_hxc input
    for i in 1:Nrmesh
        V_hxc[i,ispin] = Vpot[i,ispin] + Zed/grid.r[i]
    end

    # Solve for all states
    thresh0 = 1.0e-10
    nstop = 0
    mode = 1 # for lschps
    #XXX this is used in ascheq, using negative of Zval
    ze2 = -Zval # should be 2*Zval in Ry unit
    diff_V = Inf

    for iterSCF in 1:NiterMax

        #println("\niterSCF = ", iterSCF)

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
            if rel == 1
                Enl[iwf], nstop = lschps!( mode, Zval, thresh0, 
                    grid, nn[iwf], ll[iwf], Enl[iwf], Vpot, psi1
                )
            elseif rel == 0
                Enl[iwf], nstop = ascheq!(
                    nn[iwf], ll[iwf], Enl[iwf], grid, Vpot, ze2, thresh0, psi1, nstop
                )
            else
                error("rel=$rel is not supported")
            end
        end

        #
        # calculate charge density (spherical approximation)
        #
        fill!(Rhoe, 0.0)
        for iwf in 1:Nwf, ir in 1:Nrmesh
            # this is for ispin=1
            Rhoe[ir,ispin] += oc[iwf] * psi[ir,iwf]^2
        end
        #integRho = PWDFT.integ_simpson(Nrmesh, Rhoe[:,1], grid.rab) 
        #println("integRho = ", integRho)

        radial_poisson_solve!(0, 2, grid, Rhoe, V_h)
        #println("sum V_h = ", sum(V_h))

        Rhoe_radial[:] .= Rhoe[:] ./ grid.r2[:] ./ (4π) # 

        calc_epsxc_Vxc_VWN!(
            xc_calc, Rhoe_radial,
            epsxc,
            V_xc
        )
        #println("sum Rhoe_radial = ", sum(Rhoe_radial))
        #println("sum epsxc = ", sum(epsxc))
        #println("sum V_xc = ", sum(V_xc))

        for i in 1:Nrmesh
            V_hxc_new[i,ispin] = V_h[i] + V_xc[i,ispin]
            V_new[i,ispin] = -Zed/grid.r[i] + vxt[i] + V_hxc_new[i,ispin]
        end

        diff_V = LinearAlgebra.norm(V_hxc .- V_hxc_new)
        @printf("iterSCF = %3d diff_V = %18.10e\n", iterSCF, diff_V)
        if diff_V < 1e-10
            println("!!!! CONVERGED !!!")
            break
        end

        for i in 1:Nrmesh
            # Mix
            V_hxc[i,ispin] = 0.5*V_hxc_new[i,ispin] + 0.5*V_hxc[i,ispin]
            # Set new potential
            Vpot[i,ispin] = -Zed/grid.r[i] + V_hxc[i,ispin]
        end

    end

    println("\nAt the end of SCF")
    println("Final energy levels:")
    for iwf in 1:Nwf
        @printf("%3d %18.10f\n", iwf, Enl[iwf])
    end

    excgga = zeros(Float64, Nrmesh)
    energies = calc_energies_ae_atom(
           Zval, grid, Rhoe, V_h, V_xc, epsxc, excgga, Nwf,
           Nspin, Enl, oc
    )
    @show energies

    return

end
