push!(LOAD_PATH, pwd())

using Printf
using LAPWDFT

function create_lattice_vars()
    LatVecs = zeros(3,3)
    A = 5.13
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]
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
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    readspecies!(1, "DATA_species/Si.in", atsp_vars, mt_vars, apwlo_vars)
    readspecies!(2, "DATA_species/Pt.in", atsp_vars, mt_vars, apwlo_vars)

    init_zero!( mt_vars )

    println("before nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    checkmt!( latt_vars, atm_vars, atsp_vars.spsymb, mt_vars )
    genrmesh!( atm_vars, atsp_vars, mt_vars )
    init_packed_mtr!(mt_vars)

    println("after nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    for isp in 1:Nspecies
       
        println()
        println("species: ", isp)
        println()

        println("rsp = ")
        println(atsp_vars.rsp[isp][1:5])

        println("wmrt = ")
        println(mt_vars.wrmt[isp][1:5])
   
        println("wcmrt = ")
        println(mt_vars.wrcmt[isp][1:5])

        #println("wprmt 1 = ")
        #println(mt_vars.wprmt[1,1:5,is])

        println("wprcmt 1 = ")
        println(mt_vars.wprcmt[isp][1,1:5])
        println("wprcmt 2 = ")
        println(mt_vars.wprcmt[isp][2,1:5])
        println("wprcmt 3 = ")
        println(mt_vars.wprcmt[isp][3,1:5])
        println("wprcmt 4 = ")
        println(mt_vars.wprcmt[isp][4,1:5])

        #println("wprmt 2 = ")
        #println(mt_vars.wprmt[2,1:5,1])

        #println("wprmt 3 = ")
        #println(mt_vars.wprmt[3,1:5,1])

        #println("wprmt 4 = ")
        #println(mt_vars.wprmt[4,1:5,1])
    end

    println("lmmaxi   = ", mt_vars.lmmaxi)
    println("npmtmax  = ", mt_vars.npmtmax)
    println("npmt     = ", mt_vars.npmt)
    println("npcmtmax = ", mt_vars.npcmtmax)
    println("npcmt    = ", mt_vars.npcmt)

end

@time main()