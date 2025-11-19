# Stuffs related to states, spins, magnetic moments go here

mutable struct ElectronicChargesStates
    spinpol::Bool
    nspinor::Int64
    nempty0::Float64 #!!! This is Float64
    nempty::Int64
    #
    chgzn::Float64
    chgval::Float64
    chgcr::Vector{Float64}
    chgcrlk::Vector{Float64}
    chgcrtot::Float64
    chgexs::Float64
    chgtot::Float64
    #
    nstfv::Int64
    nstsv::Int64
    swidth::Float64
    occsv::Matrix{Float64}
    evalsv::Matrix{Float64}
    efermi::Float64
    occmax::Float64
    #
    tevecsv::Bool
end

function ElectronicChargesStates(
    atoms::Atoms, atsp_vars, Nkpt::Int64;
    nempty0=4.0, spinpol=false, chgexs=0.0,
    swidth=0.001,
    nempty::Union{Int64,Nothing}=nothing
)
    # nempty0 is number of empty states per atom. It is a Float64
    # nspinor should be read from something else?

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Nspecies = atoms.Nspecies

    spzn = atsp_vars.spzn
    occsp = atsp_vars.occsp
    nstsp = atsp_vars.nstsp
    spcore = atsp_vars.spcore

    chgzn = 0.0
    chgval = 0.0
    chgcr = zeros(Float64, Nspecies)
    chgcrlk = zeros(Float64, Natoms)
    chgcrtot = 0.0

    for ia in 1:Natoms
        isp = atm2species[ia]
        chgzn += spzn[isp]
        chgcr[isp] = 0.0
        for ist = 1:nstsp[isp]
            if spcore[isp][ist]
                chgcr[isp] += occsp[isp][ist]
            else 
                chgval += occsp[isp][ist]
            end
        end 
        chgcrtot += chgcr[isp]
    end
    # add excess charge
    chgval += chgexs
    # total charge
    chgtot = chgcrtot + chgval
  
    # nempty is integer, while nempty0 is Float64
    if isnothing(nempty)
        nempty = round(Int64, nempty0*max(Natoms,1))
        if nempty < 1
            nempty = 1
        end
    end

    if spinpol
        nspinor = 2
        occmax = 1.0
    else
        nspinor = 1
        occmax = 2.0
    end

    nstfv = round(Int64, chgval/2.0) + nempty + 1
    nstsv = nspinor*nstfv

    occsv = zeros(Float64, nstsv, Nkpt)
    evalsv = zeros(Float64, nstsv, Nkpt)
    efermi = 0.0

    if spinpol
        tevecsv = true
    else
        tevecsv = false
    end

    return ElectronicChargesStates(
        spinpol, nspinor, nempty0, nempty,
        chgzn, chgval, chgcr, chgcrlk,
        chgcrtot, chgexs, chgtot,
        nstfv, nstsv, swidth, occsv, evalsv,
        efermi, occmax,
        tevecsv
    )
end