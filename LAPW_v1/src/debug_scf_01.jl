# This is currently the entry point
# We may call @infiltrate within this function and investigate various variables

function debug_scf_01()

    # Read the input file
    elk_input = read_elk_input()
    # assign some variables
    spinpol = elk_input.spinpol
    #
    bfieldc = elk_input.bfieldc
    bfieldc0 = elk_input.bfieldc
    # NOTE: in Elk, the input variable name is bfieldc,
    #       but this is assigned to bfieldc0
    #bfcmt = elk_input.bfcmt
    #bfcmt0 = elk_input.bfcmt
    #
    spinorb = elk_input.spinorb
    spinsprl = elk_input.spinsprl
    nosource = elk_input.nosource
    cmagz = elk_input.cmagz

    ncmag = false # read from elk.in ?

    # Initialize Atoms
    atoms = create_atoms_from_elk_input(elk_input)

    # FIXME: bfcmt shoule be read from elk.in
    bfcmt = zeros(Float64, 3, atoms.Natoms)
    bfcmt0 = zeros(Float64, 3, atoms.Natoms)

    # Setup symmetry variables
    sym_vars = SymmetryVars()
    findsymlat!(sym_vars, atoms)
    findsymcrys!(sym_vars, atoms, spinpol = spinpol, bfieldc0 = bfieldc0)
    findsymsite!(sym_vars, atoms, spinpol = spinpol, bfieldc0 = bfieldc0)

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols
    specs_info = Vector{SpeciesInfo}(undef,Nspecies)
    for isp in 1:Nspecies
        specs_info[isp] = SpeciesInfo(spsymb[isp]*".in")
    end

    checkmt!(atoms, specs_info)
    # make the muffin-tin mesh commensurate with lradstp
    lradstp = 4 # in muffin tin
    for isp in 1:Nspecies
        nrmt = specs_info[isp].nrmt
        specs_info[isp].nrmt -= (nrmt-1)%lradstp
        #nrcmt[is] = (nrmt[is] - 1)/lradstp + 1
    end

    atsp_vars = AtomicSpeciesVars(atoms, specs_info)
    mt_vars = MuffinTins(specs_info, atsp_vars)
    apwlo_vars = APWLOVars(atoms, specs_info, mt_vars)

    allatoms!(atsp_vars)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Default value for rgkmax and gmaxvr
    rgkmax = 7.0
    gmaxvr = 12.0

    # use average muffin-tin radius (default)
    rsum = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        rsum += mt_vars.rmt[isp]
    end
    rsum = rsum/Natoms
    gkmax = rgkmax/rsum

    if gmaxvr <= 2.0*gkmax
        # gmaxvr is too small, set it to the minimum value
        @info "gengvec: gmaxvr will be set to 2*gkmax"
        gmaxvr = 2*gkmax
    end

    println("gmaxvr = ", gmaxvr)
    println("gkmax = ", gkmax)
    # Compute ecutrho from gmaxvr
    ecutrho = 0.5*gmaxvr^2
    # Compute ecutwfc from gkmax
    ecutwfc = 0.5*gkmax^2
    # This is required in PWGrid constructor
    dual = ecutrho/ecutwfc
    sym_info = SymmetryInfo(atoms)
    @info "sym_info.Nsyms = $(sym_info.Nsyms)"
    @info "sym_info.Nrots = $(sym_info.Nrots)"

    # FIXME: need to pass k-points information from elk_input
    @info "KPoints read from pwdftjl_kpoints.jldat"
    pw = PWGrid(
        ecutwfc, atoms.LatVecs, dual=dual,
        Ns_ = deserialize("pwdftjl_Ns.jldat"), # XXX
        #kpoints=KPoints(atoms, elk_input.ngridk, [0,0,0], sym_info.s)
        kpoints = deserialize("pwdftjl_kpoints.jldat")
    )
    println(pw)

    # check for collinearity in the z-direction and set the dimension of the
    # magnetization and exchange-correlation vector fields
    epslat  = 1e-6
    if spinpol
        ndmag = 1
        if ( abs(bfieldc0[1]) > epslat) || (abs(bfieldc0[2]) > epslat)
            ndmag = 3
            @info "ndmag is set to 3"
        end
        for ia in 1:Natoms
            if (abs(bfcmt0[1,ia]) > epslat || (abs(bfcmt0[2,ia]) > epslat))
                ndmag = 3
                @info "ndmag is set to 3"
            end
        end
        # spin-orbit coupling is non-collinear in general
        if spinorb
            ndmag = 3
        end
        # source-free fields and spin-spirals must be non-collinear
        if nosource || spinsprl
            ndmag = 3
            cmagz = false
        end
        # force collinear magnetism along the z-axis if required
        if cmagz
            ndmag = 1
        end
    else
        ndmag = 0
    end


    # Create a version for full GVectors, where Ng=Npoints
    gvec_full = GVectorsFull(pw.Ns, pw.RecVecs)
    # Some quantities for interstitial density and potentials can use gvec_full
    ffacg = genffacgp(pw, mt_vars.rmt, gvec_full = gvec_full)
    cfunig, cfunir = gencfun(pw, atoms, ffacg, gvec_full = gvec_full)
    #vsig = zeros(ComplexF64, gvec_full.Ng) # also use gvec_full ?

    #ffacg = genffacgp(pw, mt_vars.rmt, gvec_full = nothing)
    #cfunig, cfunir = gencfun(pw, atoms, ffacg, gvec_full = nothing)
    vsig = zeros(ComplexF64, pw.gvec.Ng)


    # XXX: use simpler name?
    elec_chgst = ElectronicChargesStates(
        atoms, atsp_vars, pw.gvecw.kpoints.Nkpt,
        spinpol = spinpol,
        nempty = elk_input.nempty
    )

    # Initialize rhomt and rhoir
    npmt = mt_vars.npmt
    Npoints = prod(pw.Ns)
    rhomt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        rhomt[ia] = zeros(Float64, npmt[isp])
    end
    rhoir = zeros(Float64, Npoints)
    #
    rhoinit!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir )

    # magnetization
    if spinpol
        magmt = Vector{Vector{Float64}}(undef, Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            magmt[ia] = zeros(Float64, npmt[isp])
        end
        magir = zeros(Float64, Npoints, ndmag)
        maginit!(atoms, rhomt, rhoir, ncmag, ndmag, bfcmt, bfieldc, magmt, magir)
    else
        magmt = nothing
        magir = nothing
    end

    # Hartree potential (MT and interstitial)
    vclmt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vclmt[ia] = zeros(Float64, npmt[isp])
    end
    vclir = zeros(Float64, Npoints)

    # Exchange correlation potentials
    epsxcmt = Vector{Vector{Float64}}(undef,Natoms)
    vxcmt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        epsxcmt[ia] = zeros(Float64, npmt[isp])
        vxcmt[ia] = zeros(Float64, npmt[isp])
    end
    if spinpol
        bxcmt = Vector{Matrix{Float64}}(undef,Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bxcmt[ia] = zeros(Float64, npmt[isp], ndmag)
        end
    else
        bxcmt = nothing
    end
    #
    epsxcir = zeros(Float64, Npoints)
    vxcir = zeros(Float64, Npoints)
    if spinpol
        bxcir = zeros(Float64, Npoints, ndmag)
    else
        bxcir = nothing
    end

    vsmt = Vector{Vector{Float64}}(undef, Natoms)
    vsir = zeros(Float64, Npoints)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vsmt[ia] = zeros(Float64, npmt[isp])
    end

    if spinpol
        bsir = zeros(Float64, Npoints, ndmag)
        bsmt = Vector{Matrix{Float64}}(undef,Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bsmt[ia] = zeros(Float64, npmt[isp], ndmag)
        end
    else
        bsir = nothing
        bsmt = nothing
    end

    # Old potentials (to be mixed)
    vsmt_old = Vector{Vector{Float64}}(undef, Natoms)
    vsir_old = zeros(Float64, Npoints)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vsmt_old[ia] = zeros(Float64, npmt[isp])
    end
    if spinpol
        bsir_old = zeros(Float64, Npoints, ndmag)
        bsmt_old = Vector{Matrix{Float64}}(undef,Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bsmt_old[ia] = zeros(Float64, npmt[isp], ndmag)
        end
    else
        bsir_old = nothing
        bsmt_old = nothing
    end

    core_states = CoreStatesVars(atoms, atsp_vars, mt_vars)
    apwlo_ints = APWLOIntegrals(atoms, mt_vars, apwlo_vars)
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    # This should also depend on spins. We might want to combine Nkpt and Nspin.
    nmat = zeros(Int64, Nkpt)
    for ik in 1:Nkpt
        nmat[ik] = Ngw[ik] + apwlo_vars.nlotot
        # assert ntsfv > nmat[ik]
    end
    nmatmax = maximum(nmat) # not used?
    efermi = 0.0


    # Initial call to potks?
    potks!(
        atoms, atsp_vars, mt_vars, pw, sym_vars, 
        rhomt, rhoir,
        vclmt, vclir,
        epsxcmt, epsxcir,
        vxcmt, vxcir,
        vsmt, vsir;
        magmt = magmt, magir = magir,
        bsmt = bsmt, bsir = bsir,
        bxcir = bxcir, bxcmt = bxcmt,
        cfunir = cfunir,
        bfieldc = bfieldc, bfcmt = bfcmt,
        spinpol = spinpol, ncmag = ncmag
    )


    # Fourier transform of interstitial Kohn-Sham equation
    genvsig!(pw, vsir, cfunir, vsig)
    # XXX vsig will be different from Elk result because Elk uses more G-vectors

    E_tot_old = Inf

    for iter_scf in 1:50

        println("\n\nENTER iter_scf = ", iter_scf)
        println("-------------------------------")

        gencore!(atoms, sym_vars.eqatoms, atsp_vars, mt_vars, vsmt, core_states)
        linengy!(atoms, sym_vars.eqatoms, mt_vars, vsmt, efermi, apwlo_vars)
        genapwfr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, vsmt)
        genlofr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, vsmt)
        calc_apwlo_integrals!(atoms, mt_vars, apwlo_vars, vsmt, apwlo_ints)

        # TODO: generate Hamiltonians for each kpoints, diagonalize them
        # and store the results to files (to be read later)
        # Using local variables for Hamiltonians, eigenvectors and eigenvalues
        for ik in 1:Nkpt
            gen_eigensystem!( ik,
                atoms, pw, mt_vars, apwlo_vars,
                apwlo_ints, elec_chgst,
                nmat, cfunig, vsig;
                bsmt=bsmt, bsir=bsir, ndmag=ndmag
            )
        end

        # Update occupation numbers
        occupy!(
            pw.gvecw.kpoints, apwlo_vars, elec_chgst;
            NiterMax=1000, epsocc=1e-8
        )

        # Calculate electron densities, charges, magnetizations, moments
        rhomag!(
            atoms, atsp_vars, pw, sym_vars, mt_vars, apwlo_vars, core_states, elec_chgst, cfunir, 
            rhomt, rhoir;
            magmt=magmt, magir=magir
        )

        # Save old potentials
        vsir_old[:] = vsir[:]
        for ia in 1:Natoms
            vsmt_old[ia][:] = vsmt[ia][:]
        end
        if spinpol
            bsir_old[:] = bsir[:]
            for ia in 1:Natoms
                bsmt_old[ia][:] = bsmt[ia][:]
            end
        end

        # New potential
        potks!(
            atoms, atsp_vars, mt_vars, pw, sym_vars, 
            rhomt, rhoir,
            vclmt, vclir,
            epsxcmt, epsxcir,
            vxcmt, vxcir,
            vsmt, vsir;
            magmt = magmt, magir = magir,
            bsmt = bsmt, bsir = bsir,
            bxcir = bxcir, bxcmt = bxcmt,
            cfunir = cfunir,
            bfieldc = bfieldc, bfcmt = bfcmt,
            spinpol = spinpol, ncmag = ncmag
        )

        E_tot = calc_energy_terms!(
            atoms, atsp_vars, core_states,
            pw, mt_vars, elec_chgst, ndmag,
            cfunir,
            rhomt, rhoir,
            vsmt,
            vclmt, vclir,
            epsxcmt, epsxcir, vxcmt, vxcir,
            bsmt, bsir, magmt, magir
        )
        println("E_tot = ", E_tot)

        dv_ir = sum((vsir - vsir_old).^2)/length(vsir)
        dv_mt = 0.0
        for ia in 1:Natoms
            dv_mt += sum((vsmt[ia] - vsmt_old[ia]).^2)/length(vsmt[ia])
        end
        println("dv_ir = $dv_ir dv_mt = $dv_mt")
        if spinpol
            dmag_ir = sum((bsir - bsir_old).^2)/length(bsir)
            dmag_mt = 0.0
            for ia in 1:Natoms
                dmag_mt += sum( (bsmt[ia][:] - bsmt_old[ia][:]).^2 )/length(magmt[ia])
            end                
            println("dmag_ir = $dmag_ir dmag_mt = $dmag_mt")
        end

        ΔE = abs(E_tot - E_tot_old)
        is_converged = ΔE < 1e-6
        @printf("%4d %18.10f %18.6e\n", iter_scf, E_tot, ΔE)
        if is_converged
            println("CONVERGED in total energy")
            break
        end
        E_tot_old = E_tot

        β_mix = 0.1
        vsir[:] = β_mix*vsir[:] + (1 - β_mix)*vsir_old[:]
        for ia in 1:Natoms
            vsmt[ia][:] = β_mix*vsmt[ia][:] + (1-β_mix)*vsmt_old[ia][:]
        end
        if spinpol
            bsir[:] = β_mix*bsir[:] + (1 - β_mix)*bsir_old[:]
            for ia in 1:Natoms
                bsmt[ia][:,:] = β_mix*bsmt[ia][:] + (1 - β_mix)*bsmt_old[ia][:]
            end
        end

        # Fourier transform of interstitial Kohn-Sham equation
        genvsig!(pw, vsir, cfunir, vsig)

    end # scf

    @infiltrate
    # open REPL and investigate the variables

    return

end




