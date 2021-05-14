push!(LOAD_PATH, pwd())

using Printf
using PWDFT
using LAPWDFT

function main()

    LatVecs = zeros(3,3)
    A = 5.13
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]

    atoms = Atoms(xyz_string_frac="""
    2

    Si  0.0  0.0  0.0
    Si  0.25 0.25 0.25
    """, in_bohr=true, LatVecs=LatVecs)

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    for isp in 1:Nspecies
        readspecies!(isp, "DATA_species/"*spsymb[isp]*".in", atsp_vars, mt_vars, apwlo_vars)
    end

    init_zero!( mt_vars )

    println("before nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    checkmt!( atoms, mt_vars )
    genrmesh!( atoms, atsp_vars, mt_vars )
    init_packed_mtr!( mt_vars )

    println("after nrsp    = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

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

    println("gmaxvr = ", gmaxvr)
    println("gkmax = ", gkmax)

    ecutrho = 0.5*gmaxvr^2
    ecutwfc = 0.5*gkmax^2

    dual = ecutrho/ecutwfc
    println("dual = ", dual)

    sym_info = SymmetryInfo(atoms)
    println("sym_info.Nsyms = ", sym_info.Nsyms)
    println("sym_info.Nrots = ", sym_info.Nrots)
    println("sym_info.irt   = ", sym_info.irt)

    pw = PWGrid( ecutwfc, LatVecs, dual=dual,
        kpoints=KPoints(atoms, [2,2,2], [0,0,0], sym_info.s)
    )
    println(pw)

    rhoinit!( atoms, atsp_vars, mt_vars, pw )
end

main()
