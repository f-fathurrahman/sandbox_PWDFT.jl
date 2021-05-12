push!(LOAD_PATH, pwd())

using LAPWDFT

function create_lattice_vars()
    LatVecs = zeros(3,3)
    LatVecs[1,:] = [10.0, 10.0, 0.0]
    LatVecs[2,:] = [10.0, 0.0, 10.0]
    LatVecs[3,:] = [0.0, 10.0, 10.0]
    lattice_vars = LatticeVars( LatVecs )
    return lattice_vars
end

function create_atomic_vars(lattice_vars)

    maxatoms = 200
    maxspecies = 8
    atposl = zeros(3,maxatoms,maxspecies)
    
    Nspecies = 2
    Natoms = [1,1]

    # species 1, atom 1
    atposl[:,1,1] = [0.0, 0.0, 0.0]
    # species 2, atom 1
    atposl[:,1,2] = [0.25, 0.25, 0.25]

    atomic_vars = AtomicVars(Nspecies, Natoms, atposl, lattice_vars)
end

function main()

    latt_vars = create_lattice_vars()
    atm_vars = create_atomic_vars(latt_vars)

    Nspecies = 2
    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    #apwlo_vars = APWLOVars(Nspecies, mt_vars.maxlapw) # XXX use lmaxapw instead of maxlapw?
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    readspecies!(1, "DATA_species/Si.in", atsp_vars, mt_vars, apwlo_vars)
    readspecies!(2, "DATA_species/Pt.in", atsp_vars, mt_vars, apwlo_vars)

    dmin, is, js = mtdmin( latt_vars, atm_vars, mt_vars.rmt )
    
    println("is = ", is)
    println("js = ", js)
    println("dmin = ", dmin)

    println("Pass here")

end

main()