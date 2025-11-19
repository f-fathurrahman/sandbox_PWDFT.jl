mutable struct CoreStatesVars
    spincore::Bool  # spincore is .true. if the core is to be treated as spin-polarised
    nspncr::Int64  # number of core spin-channels
    occcr::Vector{Vector{Float64}} # occupancies for core states
    evalcr::Vector{Vector{Float64}} # eigenvalues for core states
    rhocr::Vector{Matrix{Float64}} # radial charge density for core states
    rwfcr::Vector{Array{Float64,3}} # radial wavefunctions for core states
end

function CoreStatesVars(atoms, atsp_vars, mt_vars; spincore=false)
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nrmt = mt_vars.nrmt
    nstsp = atsp_vars.nstsp
    nrsp = atsp_vars.nrsp
    #
    if spincore
        nspncr = 2
    else
        nspncr = 1
    end
    #
    occcr = Vector{Vector{Float64}}(undef,Natoms)
    evalcr = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        occcr[ia] = zeros(Float64, nstsp[isp])
        evalcr[ia] = zeros(Float64, nstsp[isp])
    end
    # Copy some data from atsp_vars
    for ia in 1:Natoms
        isp = atm2species[ia]
        occcr[ia][:] = atsp_vars.occsp[isp][:]
        evalcr[ia][:] = atsp_vars.evalsp[isp][:]
    end
    #
    rhocr = Vector{Matrix{Float64}}(undef,Natoms)
    for ispin_core in 1:nspncr, ia in 1:Natoms
        isp = atm2species[ia]
        rhocr[ia] = zeros(Float64, nrmt[isp], ispin_core)
    end
    #
    rwfcr = Vector{Array{Float64,3}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        rwfcr[ia] = zeros(Float64, nrsp[isp], 2, nstsp[isp])
    end
    #
    return CoreStatesVars( spincore, nspncr, occcr, evalcr, rhocr, rwfcr )
end
