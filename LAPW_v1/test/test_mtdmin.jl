push!(LOAD_PATH, pwd())

using PWDFT
using LAPWDFT

function main()

    LatVecs = zeros(3,3)
    A = 5.0
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]

    atoms = Atoms(xyz_string_frac="""
    2

    Si  0.0   0.0   0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=LatVecs)

    Nspecies = 2
    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    readspecies!(1, "DATA_species/Si.in", atsp_vars, mt_vars, apwlo_vars)
    readspecies!(2, "DATA_species/Pt.in", atsp_vars, mt_vars, apwlo_vars)

    dmin, is, js = mtdmin(atoms, mt_vars.rmt)
    
    println("is = ", is)
    println("js = ", js)
    println("dmin = ", dmin)

    println("Pass here")

end

main()