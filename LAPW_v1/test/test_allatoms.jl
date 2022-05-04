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
    Pt  0.25 0.25 0.25
    """, in_bohr=true, LatVecs=LatVecs)

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    for isp in 1:Nspecies
        readspecies!(isp, "DATA_species/"*spsymb[isp]*".in", atsp_vars, mt_vars, apwlo_vars)
    end

    init_zero!(mt_vars)

    checkmt!(atoms, mt_vars )
    genrmesh!(atoms, atsp_vars, mt_vars)
    init_packed_mtr!( mt_vars )

    allatoms!(atsp_vars)

end

@time main()
@time main()