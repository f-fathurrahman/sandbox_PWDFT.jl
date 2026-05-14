# This is currently the entry point
# We may call @infiltrate within this function and investigate various variables

function debug_main()

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
    #bfieldc = elk_input.bfieldc
    bfieldc0 = elk_input.bfieldc
    # NOTE: in Elk, the input variable name is bfieldc,
    #       but this is assigned to bfieldc0
    #bfcmt0 = elk_input.bfcmt0
    #bfcmt = copy(bfcmt0)

    #ncmag = elk_input.ncmag
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
    specs_info = Vector{SpeciesInfo}(undef,Nspecies)
    for isp in 1:Nspecies
        specs_info[isp] = SpeciesInfo(elk_input.species_files[isp])
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
    apwlo_vars = APWLOVars(
        atoms, specs_info, mt_vars,
        nxoapwlo = elk_input.nxoapwlo
    )
    #println("nxoapwlo = ", elk_input.nxoapwlo)
    #@infiltrate
    #return

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

    ene_terms = EnergyTerms()

    E_tot = ene_terms.E_tot # should be a reference?

    gencore!(atoms, sym_vars.eqatoms, atsp_vars, mt_vars, potentials.vsmt, core_states)
    linengy!(atoms, sym_vars.eqatoms, mt_vars, potentials.vsmt, elec_chgst.efermi, apwlo_vars)
    #if iter_scf == 1
    #    break # debug
    #end
    genapwfr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, potentials.vsmt)
    genlofr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, potentials.vsmt)
    calc_apwlo_integrals!(atoms, mt_vars, apwlo_vars, potentials.vsmt, apwlo_ints)

    # TODO: generate Hamiltonians for each kpoints, diagonalize them
    # and store the results to files (to be read later)
    # Using local variables for Hamiltonians, eigenvectors and eigenvalues
    for ik in 1:Nkpt
        gen_eigensystem!( ik,
            atoms, atsp_vars, pw, mt_vars, apwlo_vars,
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

    print_info(ene_terms, prefix_str = "Final")

    #@infiltrate
    @exfiltrate()
    # open REPL and investigate the variables

    return

end



