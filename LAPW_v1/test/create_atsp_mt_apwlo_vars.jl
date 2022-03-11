function create_atsp_mt_apwlo_vars(atoms)

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
    init_packed_mtr!( mt_vars )
    #
    allatoms!(atsp_vars)

    return atsp_vars, mt_vars, apwlo_vars
end