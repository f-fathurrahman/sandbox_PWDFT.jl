include("r3frac.jl")
include("r3mv.jl")
include("LatticeVars.jl")
include("AtomicVars.jl")
include("AtomicSpeciesVars.jl")
include("MuffinTins.jl")
include("APWLOVars.jl")
include("readspecies.jl")

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
