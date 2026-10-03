
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

    V_ps_loc = zeros(Float64, Nrmesh)
    #
    # define the v_out function
    for ir in 1:ir_loc
        V_ps_loc[ir] = xc[1]*j1[ir,1] + xc[2]*j1[ir,2]
    end
    
    # Beyond rcloc V_ps_loc should be the same as Vpot
    for ir in (ir_loc+1):Nrmesh
        V_ps_loc[ir] = Vpot[ir]
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

    println("Setting appropriate energies Enls")
    nstoae = ld1x_input.nstoae
    Enls = ld1x_input.Enls
    Nwfs = ld1x_input.Nwfs
    #
    SMALL_ENERGY = 1e-13
    for iwfs in 1:Nwfs
        iwf = nstoae[iwfs]
        if abs(Enls[iwfs]) <= SMALL_ENERGY
            Enls[iwfs] = Enl[iwf]  # just assign AE energy
        end
    end

    pseudotype = 3 #XXX HARCODED
    rcut = ld1x_input.rcut
    rcutus = ld1x_input.rcutus
    Nbeta = ld1x_input.Nbeta
    fit_to_arbitrary_energy = ld1x_input.fit_to_arbitrary_energy
    #
    idx_rcut = zeros(Int64, Nbeta)
    idx_rcutus = zeros(Int64, Nbeta)
    idx_rcloc = 0
    idx_rbeta = zeros(Int64, Nbeta)
    # find ik=idx_rcut, ikus=idx_rcutus, and ikloc=idx_rcloc
    for ibeta in 1:Nbeta
        for ir in 1:Nrmesh
            if grid.r[ir] < rcut[ibeta]
                idx_rcut[ibeta] = ir
            end
            #
            if grid.r[ir] < rcutus[ibeta]
                idx_rcutus[ibeta] = ir
            end
            #
            if grid.r[ir] < rcloc
                idx_rcloc = ir
            end
        end
        # make idx_rcut odd
        if idx_rcut[ibeta]%2 == 0
            idx_rcut[ibeta] += 1
        end
        # make idx_rcutus odd
        if idx_rcutus[ibeta]%2 == 0
            idx_rcutus[ibeta] += 1
        end
        #
        if pseudotype == 3
            idx_rbeta[ibeta] = max( idx_rcutus[ibeta] + 10, idx_rcloc + 5 )
        else
            idx_rbeta[ibeta] = max( idx_rcut[ibeta] + 10, idx_rcloc + 5)
        end
        #
        @assert idx_rbeta[ibeta] < Nrmesh
    end

    lls = ld1x_input.lls
    for ibeta in 1:Nbeta, jbeta in 1:Nbeta
        if (lls[ibeta] == lls[jbeta]) && (idx_rbeta[jbeta] > idx_rbeta[ibeta])
            idx_rbeta[ibeta] = idx_rbeta[jbeta] # choose the larger
        end
    end
    # XXX: is this needed?
    irc = maximum(idx_rbeta) + 8 # offset (arbitrarily?) by 8
    if irc%2 == 0
        irc += 1
    end

    psipaw = zeros(Float64, Nrmesh, Nbeta)
    gi = zeros(Float64, Nrmesh)
    for ibeta in 1:Nbeta
        nwf0 = nstoae[ibeta]
        if fit_to_arbitrary_energy[ibeta]
            @views _set_psi_in!( Zval, Vpot, grid, idx_rcutus[ibeta], lls[ibeta], Enls[ibeta], psipaw[:,ibeta] )
            #FIXME: jjs is not yet used
        else
            #write(*,'(1x,A,I4,A)') 'ns = ', ns, ' new(ns) is false'
            ℓ = lls[ibeta]
            nst = (ℓ + 1)*2
            psipaw[:,ibeta] = psi[:,nwf0]
            #
            for ir in 1:Nrmesh
                gi[ir] = psipaw[ir,ibeta] * psipaw[ir,ibeta]
            end
            #
            nrm1 = sqrt( integ_0_inf_dr(gi, grid, Nrmesh, nst) )
            # normalize
            psipaw[:,ibeta] = psipaw[:,ibeta]/nrm1
        end
    end

    ecutrho = 0.0
    ecutwfc = 0.0
    ocs = ld1x_input.ocs
    els = ld1x_input.els
    psi_in = zeros(Float64, Nrmesh)
    psipsus = zeros(Float64, Nrmesh, Nbeta)
    phis = zeros(Float64, Nrmesh, Nbeta)
    chis = zeros(Float64, Nrmesh, Nbeta)
    for ibeta in 1:Nbeta
        ℓ = lls[ibeta]
        nst = (ℓ + 1)*2
        nwf0 = nstoae[ibeta]
        if fit_to_arbitrary_energy[ibeta]
            occ = 1.0
        else
            occ = ocs[ibeta]
        end
        #
        # save the all-electron function for the PAW setup
        psi_in[1:Nrmesh] = psipaw[1:Nrmesh,ibeta]
        #
        # compute the phi functions
        psipsus[:,ibeta] = psi_in[:]
        if idx_rcutus[ibeta] != idx_rcut[ibeta]
            #
            @views compute_phius!(grid, ℓ, idx_rcutus[ibeta], psipsus[:,ibeta], phis[:,ibeta], xc, 1, els[ibeta])
            ecutwfc = max(ecutwfc, 2.0*xc[5]^2)
            println("ecutwfc = $ecutwfc Ry")
            lbes4 = true
        end
        println("lbes4 = ", lbes4)
        @views compute_chi!(
            grid, V_ps_loc, ℓ, idx_rbeta[ibeta],
            phis[:,ibeta], chis[:,ibeta], xc, Enls[ibeta], lbes4
        )

    end

    bmat = zeros(Float64, Nbeta, Nbeta)
    for ibeta in 1:Nbeta, jbeta in 1:Nbeta
        if lls[ibeta] == lls[jbeta] # also need check jjs && abs(jjs(ns)-jjs(ns1)) < 1.e-7_dp )
            nst = (lls[ibeta] + 1)*2
            idx_r = idx_rbeta[ibeta]
            for ir in 1:Nrmesh
                gi[ir] = phis[ir,ibeta]*chis[ir,jbeta]
            end
            bmat[ibeta,jbeta] = integ_0_inf_dr(gi, grid, idx_r, nst)
            println("ibeta=$ibeta jbeta=$jbeta bmat=$(bmat[ibeta,jbeta])")
        end
    end
    println("bmat = ")
    display(bmat); println()

    B = zeros(Float64, Nbeta, Nbeta)
    for ibeta in 1:Nbeta, jbeta in 1:Nbeta
        B[ibeta,jbeta] = bmat[ibeta,jbeta]
    end
    #
    # compute the inverse of the matrix B_{ij}:  B_{ij}^-1
    Binv = inv(B)
    #println("Binv = ")
    #display(Binv); println()

    # compute the beta functions
    beta_prj = zeros(Float64, Nrmesh, Nbeta)
    for ibeta in 1:Nbeta, jbeta in 1:Nbeta
        for ir in 1:Nrmesh
            beta_prj[ir,ibeta] += Binv[jbeta,ibeta]*chis[ir,jbeta]
        end
    end
    # B and Binv are not used anymore

    # the following is only for pseudotype == 3
    #
    # compute the Q functions
    qvan = zeros(Float64, Nrmesh, Nbeta, Nbeta)
    qq = zeros(Float64, Nbeta, Nbeta)
    for ibeta in 1:Nbeta, jbeta in 1:ibeta
        idx_rbeta_max = max(idx_rbeta[ibeta], idx_rbeta[jbeta])
        #XXX relativistic case is not covered here
        for ir in 1:idx_rbeta_max
            qvan[ir,ibeta,jbeta] = psipsus[ir,ibeta]*psipsus[ir,jbeta] - phis[ir,ibeta]*phis[ir,jbeta]
            gi[ir] = qvan[ir,ibeta,jbeta]
        end
        for ir in idx_rbeta_max+1:Nrmesh
            qvan[ir,ibeta,jbeta] = 0.0
        end
        #
        # and puts its integral in qq
        if lls[ibeta] == lls[jbeta]  # XXX also need to check for jjs in case of 
            nst = (lls[ibeta] + 1)*2
            qq[ibeta,jbeta] = integ_0_inf_dr(gi, grid, idx_rbeta_max, nst)
        end
        #
        # set the bmat with the eigenvalue part
        #
        bmat[ibeta,jbeta] += Enls[jbeta]*qq[ibeta,jbeta]*2 #XXX Convert to Ry ???
        #
        # Use symmetry of the n,ns1 indices to set qvan and qq and bmat
        if ibeta != jbeta
            for ir in 1:Nrmesh
                qvan[ir,jbeta,ibeta] = qvan[ir,ibeta,jbeta]
            end
            qq[jbeta,ibeta] = qq[ibeta,jbeta]
            bmat[jbeta,ibeta] += Enls[ibeta]*qq[jbeta,ibeta]*2 #XXX Convert to Ry ???
        end
    end
    
    println("The bmat + epsilon qq matrix ")
    display(bmat); println()
    println("qq matrix ")
    display(qq); println()

    #Is it the same for all spin components? 
    ddd = zeros(Float64, Nbeta, Nbeta, Nspin)
    for ispin in 1:Nspin
        ddd[:,:,ispin] = bmat[:,:]
    end

    lmx = 3
    lmx2 = 2*lmx # XXX HARCODED
    qvanl = OffsetArray(
        zeros(Nrmesh, Nbeta, Nbeta, lmx2+1),
        1:Nrmesh, 1:Nbeta, 1:Nbeta, 0:lmx2
    )
    calc_pseudo_q!(ld1x_input, grid, qvan, qvanl, idx_rbeta)

    Nwfts = 0
    for iwfs in 1:Nwfs
        if (ocs[iwfs]  > 0.0) || ( ocs[iwfs] == 0.0 && Enls[iwfs] == 0.0 )
            Nwfts += 1
        end
    end

    # copy states used in the PP generation to testing configuration
    # Only bound states must be copied. Note that this WILL NOT WORK
    # if bound states are not used in the generation of the PP
    elts = Vector{String}(undef, Nwfts)
    nnts = zeros(Int64, Nwfts)
    llts = zeros(Int64, Nwfts)
    octs = zeros(Float64, Nwfts)
    iswts = zeros(Int64, Nwfts)
    #jjts = zeros(Float64, Nwfts)
    iwfts = 0
    nns = ld1x_input.nns
    isws = ld1x_input.isws
    #jjs = ld1x_input.jjs
    for iwfs in 1:Nwfs
        if (ocs[iwfs]  > 0.0) || ( ocs[iwfs] == 0.0 && Enls[iwfs] == 0.0 )
            iwfts += 1
            elts[iwfts] = els[iwfs]
            nnts[iwfts] = nns[iwfs]
            llts[iwfts] = lls[iwfs]
            octs[iwfts] = ocs[iwfs]
            iswts[iwfts]= isws[iwfs]
            #jjts[iwfts]= jjs[iwfs]
        end
    end

    nstoaets = zeros(Int64, Nwfts)
    el = ld1x_input.el
    for iwfts in 1:Nwfts
        nstoaets[iwfts] = 0
        for iwf in 1:Nwf
            # XXX Spin-unpolarized
            if elts[iwfts] == el[iwf]  # need to check for jjs
                nstoaets[iwfts] = iwf
            end
        end
    end

    Enlts = zeros(Float64, Nwfts)
    for iwfts in 1:Nwfts
        Enlts[iwfts] = Enl[nstoaets[iwfts]]
    end

    phits = zeros(Float64, Nrmesh, Nwfts)
    # ddd should not be modified here
    # what about qq? beta_prj ?
    #println("ddd before = "); display(ddd); println()
    for iwfts in 1:Nwfts
        #XXX convert energy and potential to Ry
        @views Enlts[iwfts] = ascheqps_Ry!(
            nnts[iwfts], llts[iwfts], 2*Enlts[iwfts], grid, 2*V_ps_loc, phits[:,iwfts],
            beta_prj, ddd, qq, lls, idx_rbeta
        )
        Enlts[iwfts] *= 0.5 # scale back to Ha
        println("Output energy (in Ha) = ", Enlts[iwfts])
        #
        @views l1dx_normalize!(ld1x_input, grid, idx_rcut, qq, beta_prj, phits[:,iwfts], llts[iwfts])
        #
    end
    
    println("\nComputing bmat")
    println("bmat before = "); display(bmat); println()
    # this array should be local
    vaux = zeros(Float64, Nrmesh, 2) # 2 for spin?
    for ibeta in 1:Nbeta, jbeta in 1:ibeta
        if lls[ibeta] == lls[jbeta] # abs(jjs(ib)-jjs(jb)) < 1.e-7_dp
            ℓ = lls[ibeta]
            nst = (ℓ + 1)*2
            #
            idx_r = idx_rbeta[ibeta]
            # This is for PAW?
            #if which_augfun == :PSQ
            for ir = 1:idx_rcut[ibeta]
                vaux[ir,1] = qvanl[ir,ibeta,jbeta,0]*V_ps_loc[ir]*2 # to Ry
            end
            #ELSE
            #for ir in 1:idx_r
            #    vaux[ir,1] = qvan[ir,ibeta,jbeta]*V_ps_loc[ir] # to Ry?
            #end
            println()
            println("ibeta=$ibeta, jbeta=$jbeta")
            println("idx_r = ", idx_r)
            println("sum qvan[1:idx_r,ibeta,jbeta] = ", sum(qvan[1:idx_r,ibeta,jbeta]))
            println("sum qvanl[1:idx_r,ibeta,jbeta,0] = ", sum(qvanl[1:idx_r,ibeta,jbeta,0]))
            println("sum V_ps_loc[1:idx_r] in Ry = ", 2*sum(V_ps_loc[1:idx_r]))
            println("sum vaux[1:idx_r,1] = ", sum(vaux[1:idx_r,1]))
            println("integral result = ", integ_0_inf_dr(vaux[1:idx_r,1], grid, idx_r, nst))
            # convert to Ry?
            bmat[ibeta,jbeta] -= integ_0_inf_dr(vaux[:,1], grid, idx_r, nst)
        end
        bmat[jbeta,ibeta] = bmat[ibeta,jbeta]
    end
    println("The ddd matrix (bmat) after descreening D coefs")
    display(bmat); println()

    rhos = zeros(Float64, Nrmesh, Nspin)
    calc_charge_ps!(
        grid, rhos, phits, Nwfts, llts, octs,
        beta_prj, Nbeta, lls, idx_rbeta, qvan, qvanl;
        Nspin = Nspin
    )

    # This will calculate new potential
    radial_poisson_solve!(0, 2, grid, rhos, V_h)
    #
    #Rhoe_radial[:] .= rhos[:] ./ grid.r2[:] ./ (4π) # using rhos
    @. Rhoe_radial[:] = (rhos + rhoc) / grid.r2 / (4π) # using rhos+rhoc
    #XXX also need rho core
    calc_epsxc_Vxc_VWN!(
        xc_calc, Rhoe_radial,
        epsxc,
        Vxc
    )
    ispin = 1
    for i in 1:Nrmesh
        vaux[i,ispin] = V_h[i] + Vxc[i,ispin]
    end

    V_ps_tot = zeros(Float64, Nrmesh, Nspin)
    for ir in 1:Nrmesh
        V_ps_tot[ir,1] = V_ps_loc[ir]
        V_ps_loc[ir] -= vaux[ir,1] # subtract vaux from V_ps_loc
    end


    @infiltrate

    return
end
