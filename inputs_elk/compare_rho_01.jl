using OffsetArrays
using Statistics: mean
using Serialization

const ELK_DATA_DIR = "./TEMP_datadir"

function get_elk_rho(atoms, mt_vars)

    raw_rhomt = deserialize(joinpath(ELK_DATA_DIR, "rhomt.jldat"))
    raw_rhoir = deserialize(joinpath(ELK_DATA_DIR, "rhoir.jldat"))
    Natoms = atoms.Natoms
    npmt = mt_vars.npmt
    elk_rhomt = Vector{Vector{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        isp = atoms.atm2species[ia]
        elk_rhomt[ia] = zeros(Float64, npmt[isp])
        elk_rhomt[ia][:] = raw_rhomt[1:npmt[isp],ia] # use ias ???
    end
    elk_rhoir = copy(raw_rhoir) # no need to reshape?
    return elk_rhomt, elk_rhoir
end

