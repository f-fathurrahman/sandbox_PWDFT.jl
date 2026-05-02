

function debug_scf_02(; NiterSCFMax = 100)

    # Prepare for temporary directory for saving some data
    tmpdir = "./tmp"
    @info "Using $tmpdir for temporary directory"
    if !isdir(tmpdir)
        mkdir(tmpdir)
        @info "$tmpdir is not yet exist. Creating it."
    end

    # Read the input file
    elk_input = read_elk_input()
    # assign some variables
    spinpol = elk_input.spinpol
    #
    bfieldc = elk_input.bfieldc
    bfieldc0 = elk_input.bfieldc
    # NOTE: in Elk, the input variable name is bfieldc,
    #       but this is assigned to bfieldc0
    bfcmt0 = elk_input.bfcmt0
    bfcmt = copy(bfcmt0)

    ncmag = elk_input.ncmag
    ndmag = elk_input.ndmag

    rgkmax = elk_input.rgkmax
    gmaxvr = elk_input.gmaxvr

    # Initialize Atoms
    atoms = create_atoms_from_elk_input(elk_input)

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
    lradstp = elk_input.lradstp
    for isp in 1:Nspecies
        nrmt = specs_info[isp].nrmt
        specs_info[isp].nrmt -= (nrmt-1)%lradstp
    end

    atsp_vars = AtomicSpeciesVars(atoms, specs_info)
    mt_vars = MuffinTins(specs_info, atsp_vars, lradstp = lradstp)
    apwlo_vars = APWLOVars(atoms, specs_info, mt_vars)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
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

    # Compute ecutrho from gmaxvr
    ecutrho = 0.5*gmaxvr^2
    # Compute ecutwfc from gkmax
    ecutwfc = 0.5*gkmax^2
    # This is required in PWGrid constructor
    dual = ecutrho/ecutwfc
    sym_info = SymmetryInfo(atoms)

    # FIXME: need to pass k-points information from elk_input
    if isfile("./pwdftjl_kpoints.jldat") && isfile("./pwdftjl_kpoints.jldat")
        @info "KPoints read from pwdftjl_kpoints.jldat"
        @info "Ns is read from pwdftjl_Ns.jldat"
        pw = PWGrid(
            ecutwfc, atoms.LatVecs, dual=dual,
            Ns_ = deserialize("pwdftjl_Ns.jldat"), # XXX
            kpoints = deserialize("pwdftjl_kpoints.jldat")
        )
    else
        pw = PWGrid(
            ecutwfc, atoms.LatVecs, dual=dual,
            kpoints = KPoints(atoms, elk_input.ngridk, [0,0,0], sym_info.s)
        )
    end
    println(pw)


    need_gvec_full = true
    if need_gvec_full
        # Create a version for full GVectors, where Ng=Npoints
        gvec_full = GVectorsFull(pw.Ns, pw.RecVecs)
        # Some quantities for interstitial density and potentials can use gvec_full
        ffacg = genffacgp(pw, mt_vars.rmt, gvec_full = gvec_full)
        cfunig, cfunir = gencfun(pw, atoms, ffacg, gvec_full = gvec_full)
        #vsig = zeros(ComplexF64, gvec_full.Ng) # also use gvec_full ?
    else
        ffacg = genffacgp(pw, mt_vars.rmt, gvec_full = nothing)
        cfunig, cfunir = gencfun(pw, atoms, ffacg, gvec_full = nothing)
    end

    #XXX Call this here?
    allatoms!(atsp_vars)

    # XXX: use simpler name?
    elec_chgst = ElectronicChargesStates(
        atoms, atsp_vars, pw.gvecw.kpoints.Nkpt,
        spinpol = spinpol,
        nempty = elk_input.nempty
    )

    # Initialize rhomt and rhoir
    densities = LAPWDensities(atoms, atsp_vars, mt_vars, pw, elk_input)
    potentials = LAPWPotentials(atoms, mt_vars, pw, elk_input)


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

    # Initial call to potks?
    potks!(
        atoms, atsp_vars, mt_vars, pw, sym_vars, elk_input, 
        cfunir, densities, potentials
    )

    # Fourier transform of interstitial Kohn-Sham equation
    genvsig!(pw, potentials.vsir, cfunir, potentials.vsig)
    # XXX vsig will be different from Elk result because Elk uses more G-vectors


    Npoints = prod(pw.Ns)
    npmt = mt_vars.npmt
    # Old potentials (to be mixed)
    vsmt_old = Vector{Vector{Float64}}(undef, Natoms)
    vsir_old = zeros(Float64, Npoints)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vsmt_old[ia] = zeros(Float64, npmt[isp])
    end
    if spinpol
        bsir_old = zeros(Float64, Npoints, ndmag)
        bsmt_old = Vector{Matrix{Float64}}(undef, Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bsmt_old[ia] = zeros(Float64, npmt[isp], ndmag)
        end
    else
        bsir_old = nothing
        bsmt_old = nothing
    end

    ene_terms = EnergyTerms()
    ene_terms_old = EnergyTerms()

    E_tot = ene_terms.E_tot # should be a reference?
    E_tot_old = Inf
    etot_conv_thr = 1e-4
    dv_conv_thr = 1e-6
    Nconv = 0

    betamix = 0.1
    mixer = BroydenMixer_LAPW(potentials.vsmt, potentials.vsir, betamix, mixdim = 4)
    if spinpol
        mixer_b = BroydenMixer_LAPW(potentials.bsmt, potentials.bsir, betamix, mixdim = 4)
    end

    for iter_scf in 1:NiterSCFMax

        println("\n\nENTER iter_scf = ", iter_scf)
        println("-------------------------------")

        gencore!(atoms, sym_vars.eqatoms, atsp_vars, mt_vars, potentials.vsmt, core_states)
        linengy!(atoms, sym_vars.eqatoms, mt_vars, potentials.vsmt, elec_chgst.efermi, apwlo_vars)
        genapwfr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, potentials.vsmt)
        genlofr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, potentials.vsmt)
        calc_apwlo_integrals!(atoms, mt_vars, apwlo_vars, potentials.vsmt, apwlo_ints)

        # TODO: generate Hamiltonians for each kpoints, diagonalize them
        # and store the results to files (to be read later)
        # Using local variables for Hamiltonians, eigenvectors and eigenvalues
        for ik in 1:Nkpt
            gen_eigensystem!( ik,
                atoms, pw, mt_vars, apwlo_vars,
                apwlo_ints, elec_chgst,
                nmat, cfunig, potentials.vsig;
                bsmt = potentials.bsmt, bsir = potentials.bsir, ndmag = ndmag
            )
        end

        # Update occupation numbers
        occupy!(
            pw.gvecw.kpoints, apwlo_vars, elec_chgst;
            NiterMax = 1000, epsocc = 1e-8
        )

        # Calculate electron densities, charges, magnetizations, moments
        rhomag!(
            atoms, atsp_vars, pw, sym_vars, mt_vars, apwlo_vars, core_states, elec_chgst, cfunir, 
            densities.rhomt, densities.rhoir;
            magmt = densities.magmt, magir = densities.magir
        )

        # Save old potentials
        vsir_old[:] = potentials.vsir[:]
        for ia in 1:Natoms
            vsmt_old[ia][:] = potentials.vsmt[ia][:]
        end
        if spinpol
            bsir_old[:] = potentials.bsir[:]
            for ia in 1:Natoms
                bsmt_old[ia][:] = potentials.bsmt[ia][:]
            end
        end

        # New potential
        potks!(
            atoms, atsp_vars, mt_vars, pw, sym_vars, elk_input, 
            cfunir, densities, potentials
        )

        calc_energy_terms!(
            ene_terms,
            atoms, atsp_vars, core_states,
            pw, mt_vars, elec_chgst, ndmag,
            cfunir,
            densities.rhomt, densities.rhoir,
            potentials.vsmt,
            potentials.vclmt, potentials.vclir,
            potentials.epsxcmt, potentials.epsxcir,
            potentials.vxcmt, potentials.vxcir,
            potentials.bsmt, potentials.bsir,
            densities.magmt, densities.magir,
            potentials.bxcmt, potentials.bxcir
        )
        E_tot = ene_terms.E_tot # new calc energy
        #println("E_tot = ", E_tot)

        dv_ir = sum(abs.(potentials.vsir - vsir_old))/length(potentials.vsir)
        dv_mt = 0.0
        for ia in 1:Natoms
            dv_mt += sum(abs.(potentials.vsmt[ia] - vsmt_old[ia]))/length(potentials.vsmt[ia])
        end
        dv_mt /= Natoms # Normalize by no. of atoms
        println("dv_ir = $dv_ir dv_mt = $dv_mt")
        if spinpol
            dmag_ir = sum(abs.(potentials.bsir - bsir_old))/length(potentials.bsir)
            dmag_mt = 0.0
            for ia in 1:Natoms
                dmag_mt += sum( abs.(potentials.bsmt[ia][:] - bsmt_old[ia][:]) )/length(potentials.bsmt[ia])
            end
            dmag_mt /= Natoms # Normalize by no. of atoms
            println("dmag_ir = $dmag_ir dmag_mt = $dmag_mt")
        end
        dv = dv_mt + dv_ir
        if dv <= dv_conv_thr
            is_converged_dv = true
        else
            is_converged_dv = false
        end
        println("is_converged_dv = ", is_converged_dv)

        ΔE = abs(E_tot - E_tot_old)
        if ΔE <= etot_conv_thr
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        # 
        is_converged = (Nconv >= 2) && is_converged_dv

        #ΔE_terms = abs(ene_terms - ene_terms_old)
        #print_info(ene_terms)
        #print_info(ΔE_terms, prefix_str="Diff ")

        @printf("%4d %18.10f %18.6e\n", iter_scf, E_tot, ΔE)
        if is_converged
            println("CONVERGED in total energy (convergence achieved two times in row)")
            println("             and total potential")
            break
        end
        E_tot_old = E_tot
        ene_terms_old = copy(ene_terms)

        # This is using BroydenMixer_LAPW
        do_mix_LAPW!(mixer, potentials.vsir, vsir_old, potentials.vsmt, vsmt_old, iter_scf)
        println("Using BroydenMixer_LAPW")
        if spinpol
            do_mix_LAPW!(mixer_b, potentials.bsir, bsir_old, potentials.bsmt, bsmt_old, iter_scf)
        end

        # Fourier transform of interstitial Kohn-Sham equation
        genvsig!(pw, potentials.vsir, cfunir, potentials.vsig)

    end # scf

    print_info(ene_terms, prefix_str = "Final")

    @infiltrate
    # open REPL and investigate the variables

    return

end

