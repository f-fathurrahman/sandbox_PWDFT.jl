# This is currently the entry point
# We may call @infiltrate within this function and investigate various variables

function debug_main()

    # Read the input file
    elk_input = read_elk_input()
    
    # Initialize Atoms
    atoms = create_atoms_from_elk_input(elk_input)

    # Setup symmetry variables
    sym_vars = SymmetryVars()
    findsymlat!(sym_vars, atoms)
    findsymcrys!(sym_vars, atoms)
    findsymsite!(sym_vars, atoms)

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

    # FIXME: need to pass k-points infor from elk_input
    pw = PWGrid(
        ecutwfc, atoms.LatVecs, dual=dual,
        kpoints=KPoints(atoms, elk_input.ngridk, [0,0,0], sym_info.s)
    )
    println(pw)

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

    potxcmt!(atoms, mt_vars, rhomt, epsxcmt, vxcmt)

    epsxcir = zeros(Float64, Npoints)
    vxcir = zeros(Float64, Npoints)
    potxcir!(rhoir, epsxcir, vxcir)

    # Symmetrize
    symrfmt!(atoms, mt_vars, sym_vars, vxcmt)
    symrfir!(pw, sym_vars, vxcir)

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
    
    # generating the effective magnetic fields is skipped

    # generating the tau-DFT effective potential is skipped
    
    # .... This is the end of potks


    # Create a version for full GVectors, where Ng=Npoints

    gvec_full = GVectorsFull(pw.Ns, pw.RecVecs)

    ffacg = genffacgp(pw, mt_vars.rmt, gvec_full=gvec_full)
    cfunig, cfunir = gencfun(pw, atoms, ffacg, gvec_full=gvec_full)
    vsig = zeros(ComplexF64, pw.gvec.Ng)
    genvsig!(pw, vsir, cfunir, vsig)
    # XXX vsig will be different from Elk result because Elk uses more G-vectors
    
    core_states = CoreStatesVars(atoms, atsp_vars, mt_vars)
    gencore!(atoms, sym_vars.eqatoms, atsp_vars, mt_vars, vsmt, core_states)

    efermi = 0.0
    linengy!(atoms, sym_vars.eqatoms, mt_vars, vsmt, efermi, apwlo_vars)

    genapwfr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, vsmt)
    genlofr!(atoms, sym_vars.eqatoms, mt_vars, apwlo_vars, vsmt)


    nlomax = apwlo_vars.nlomax
    apwordmax = apwlo_vars.apwordmax
    lmaxapw = mt_vars.lmaxapw
    lmmaxo = mt_vars.lmmaxo

    oalo = Vector{Matrix{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        oalo[ia] = zeros(Float64, apwordmax, nlomax)
    end

    ololo = Vector{Matrix{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        ololo[ia] = zeros(Float64, nlomax, nlomax)
    end

    haa = Vector{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, Natoms)
    for ia in 1:Natoms
        haa[ia] = OffsetArray(
            zeros(Float64, lmmaxo, apwordmax, lmaxapw+1, apwordmax, lmaxapw+1),
            1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:apwordmax, 0:lmaxapw
        )
    end

    hloa = Vector{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, Natoms)
    for ia in 1:Natoms
        hloa[ia] = OffsetArray(
            zeros(Float64, lmmaxo, apwordmax, lmaxapw+1, nlomax),
            1:lmmaxo, 1:apwordmax, 0:lmaxapw, 1:nlomax
        )
    end

    hlolo = Vector{Array{Float64,3}}(undef, Natoms)
    for ia in 1:Natoms
        hlolo[ia] = zeros(Float64, lmmaxo, nlomax, nlomax)
    end

    @infiltrate
    # open REPL and investigate the variables

    return

end

function init_gntyry(mt_vars)
    lmmaxo = mt_vars.lmmaxo
    lmmaxapw = mt_vars.lmmaxapw
    lmaxo = mt_vars.lmaxo
    lmaxapw = mt_vars.lmaxapw
    idxlm = mt_vars.idxlm
    gntyry = zeros(ComplexF64, lmmaxo, lmmaxapw, lmmaxapw)
    for l1 in 0:lmaxapw, m1 in -l1:l1
        lm1 = idxlm[l1,m1]
        for l3 in 0:lmaxapw, m3 in -l3:l3
            lm3 = idxlm[l3,m3]
            for l2 in 0:lmaxo, m2 in -l2:l2
                lm2 = idxlm[l2,m2]
                gntyry[lm2,lm3,lm1] = gauntyry(l1,l2,l3,m1,m2,m3)
            end
        end
    end
    return gntyry
end

