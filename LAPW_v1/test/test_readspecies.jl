push!(LOAD_PATH, pwd())

using LAPWDFT

function main()
    Nspecies = 2
    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    #apwlo_vars = APWLOVars(Nspecies, mt_vars.maxlapw) # XXX use lmaxapw instead of maxlapw?
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    readspecies!(1, "DATA_species/Si.in", atsp_vars, mt_vars, apwlo_vars)
    readspecies!(2, "DATA_species/Pt.in", atsp_vars, mt_vars, apwlo_vars)

    println(atsp_vars.spsymb)
    println(atsp_vars.nstsp)
end

main()
