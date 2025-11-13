# This is currently the entry point
# We may call @infiltrate within this function and investigate various variables

function debug_main()

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
    ffacg = genffacgp(pw, mt_vars.rmt, gvec_full=gvec_full)
    cfunig, cfunir = gencfun(pw, atoms, ffacg, gvec_full=gvec_full)
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
    end


    # .... This is starting point of potks
    # XXX Need to wrap into functions, preallocate all arrays

    # Solver Hartree equation (compute electrostatic Coulomb potential)
    vclmt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vclmt[ia] = zeros(Float64, npmt[isp])
    end
    vclir = zeros(Float64, Npoints)
    potcoul!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir, vclmt, vclir )

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
    end

    if spinpol
        potxcmt!(atoms, mt_vars, rhomt, magmt, epsxcmt, vxcmt, bxcmt)
    else
        potxcmt!(atoms, mt_vars, rhomt, epsxcmt, vxcmt)
    end

    epsxcir = zeros(Float64, Npoints)
    vxcir = zeros(Float64, Npoints)
    if spinpol
        bxcir = zeros(Float64, Npoints, ndmag)
    end

    if spinpol
        potxcir!(rhoir, magir, epsxcir, vxcir, bxcir)
    else
        potxcir!(rhoir, epsxcir, vxcir)
    end

    # Symmetrize
    symrfmt!(atoms, mt_vars, sym_vars, vxcmt)
    symrfir!(pw, sym_vars, vxcir)
    if spinpol
        symrfmt!(atoms, mt_vars, sym_vars, bxcmt)
        symrfir!(pw, sym_vars, bxcir)
    end

    # effective potential from sum of Coulomb and exchange-correlation potentials
    vsmt = Vector{Vector{Float64}}(undef, Natoms)
    vsir = zeros(Float64, Npoints)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vsmt[ia] = zeros(Float64, npmt[isp])
    end
    for ia in 1:Natoms
        @views vsmt[ia][:] .= vclmt[ia][:] .+ vxcmt[ia][:]
    end
    vsir[:] = vclir[:] + vxcir[:]

    # smoothing vsir is skipped (default is zero)
    
    # Generate the effective magnetic fields
    if spinpol
        bsir = zeros(Float64, Npoints, ndmag)
        bsmt = Vector{Matrix{Float64}}(undef,Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bsmt[ia] = zeros(Float64, npmt[isp], ndmag)
        end
        genbs!( atoms, mt_vars, cfunir, ncmag,
            bfcmt, bfieldc, bxcmt, bxcir, bsmt, bsir
        )
    end


    # generating the tau-DFT effective potential is skipped
    
    # .... This is the end of potks


    # Fourier transform of interstitial Kohn-SHam equation
    genvsig!(pw, vsir, cfunir, vsig)
    # XXX vsig will be different from Elk result because Elk uses more G-vectors

    
    core_states = CoreStatesVars(atoms, atsp_vars, mt_vars)
    gencore!(atoms, sym_vars.eqatoms, atsp_vars, mt_vars, vsmt, core_states)

    efermi = 0.0
    linengy!(atoms, sym_vars.eqatoms, mt_vars, vsmt, efermi, apwlo_vars)

    genapwfr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, vsmt)
    genlofr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, vsmt)

    apwlo_ints = APWLOIntegrals( atoms, mt_vars, apwlo_vars )
    calc_apwlo_integrals!(atoms, mt_vars, apwlo_vars, vsmt, apwlo_ints)

    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    # This should also depend on spins. We might want to combine Nkpt and Nspin.
    nmat = zeros(Int64, Nkpt)
    for ik in 1:Nkpt
        nmat[ik] = Ngw[ik] + apwlo_vars.nlotot
        # assert ntsfv > nmat[ik]
    end
    nmatmax = maximum(nmat) # not used?


    # TODO: generate Hamiltonians for each kpoints, diagonalize them
    # and store the results to files (to be read later)
    # Using local variables for Hamiltonians, eigenvectors and eigenvalues
    ispin = 1
    for ik in 1:Nkpt
        gen_eigensystem( ispin, ik,
            atoms, pw, mt_vars, apwlo_vars,
            apwlo_ints, elec_chgst,
            nmat, cfunig, vsig
        )
    end

#=
    # Update occupation numbers
    occupy!(
        pw.gvecw.kpoints, apwlo_vars, elec_chgst;
        NiterMax=1000, epsocc=1e-8
    )
=#

    @infiltrate
    # open REPL and investigate the variables

    return

end




