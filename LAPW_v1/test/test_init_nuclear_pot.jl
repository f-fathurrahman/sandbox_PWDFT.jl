using Printf
using PWDFT: Atoms
using LAPWDFT

include("create_atoms.jl")

function main()

    #atoms = create_Si_fcc()
    atoms = create_SiPt_fcc()

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    for isp in 1:Nspecies
        readspecies!(isp, "DATA_species/"*spsymb[isp]*".in", atsp_vars, mt_vars, apwlo_vars)
    end

    init_zero!( mt_vars )
    checkmt!( atoms, mt_vars )
    genrmesh!( atoms, atsp_vars, mt_vars )
    init_packed_mtr!(mt_vars)

    println("Before size(atsp_vars.vcln) = ", size(atsp_vars.vcln))
    init_nuclear_pot!( atsp_vars )
    println("After  size(atsp_vars.vcln) = ", size(atsp_vars.vcln))

    for is in 1:Nspecies
        println(atsp_vars.vcln[is][1:5])
    end
end

@time main()