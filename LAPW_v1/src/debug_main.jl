function debug_main()

    elk_input = read_elk_input()
    atoms = create_atoms_from_elk_input(elk_input)

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    for isp in 1:Nspecies
        readspecies!(isp, spsymb[isp]*".in", atsp_vars, mt_vars, apwlo_vars)
    end

    init_zero!( mt_vars )
    checkmt!( atoms, mt_vars )
    genrmesh!( atoms, atsp_vars, mt_vars )
    init_packed_mtr!( mt_vars )

    init_nuclear_pot!( atsp_vars )
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
        rsum = rsum + mt_vars.rmt[isp]
    end
    rsum = rsum/Natoms
    gkmax = rgkmax/rsum

    if gmaxvr <= 2.0*gkmax
        println("INFO gengvec: gmaxvr will be set to 2*gkmax")
        gmaxvr = 2*gkmax
    end

    println("gmaxvr = ", gmaxvr)
    println("gkmax = ", gkmax)

    ecutrho = 0.5*gmaxvr^2
    ecutwfc = 0.5*gkmax^2

    dual = ecutrho/ecutwfc
    println("dual = ", dual)

    sym_info = SymmetryInfo(atoms)
    println("sym_info.Nsyms = ", sym_info.Nsyms)
    println("sym_info.Nrots = ", sym_info.Nrots)

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

    vclmt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vclmt[ia] = zeros(Float64, npmt[isp])
    end
    vclir = zeros(Float64, Npoints)
    potcoul!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir, vclmt, vclir )

    @printf("sum rhoir = %18.10f\n", sum(rhoir))
    ss = 0.0
    for ia in 1:Natoms
        ss += sum(rhomt[ia])
    end
    @printf("sum rhomt = %18.10f\n", ss)

    return
end

