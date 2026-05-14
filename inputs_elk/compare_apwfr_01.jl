using OffsetArrays
using Statistics: mean
using Serialization

const ELK_DATA_DIR = "./TEMP_datadir"


function get_elk_apwfr(atoms, mt_vars, apwlo_vars)

    raw_apwfr = deserialize(joinpath(ELK_DATA_DIR, "apwfr.jldat"))

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    lmaxapw = mt_vars.lmaxapw
    apword = apwlo_vars.apword
    APWFR_ELTYPE = OffsetVector{Vector{Matrix{Float64}}, Vector{Vector{Matrix{Float64}}}}
    elk_apwfr = Vector{APWFR_ELTYPE}(undef, Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = mt_vars.nrmt[isp]
        elk_apwfr[ia] = OffsetArray( Vector{Vector{Matrix{Float64}}}(undef, lmaxapw+1), 0:lmaxapw )
        for l in 0:lmaxapw
            elk_apwfr[ia][l] = Vector{Matrix{Float64}}(undef, apword[isp][l])
            for io in 1:apword[isp][l]
                elk_apwfr[ia][l][io] = zeros(Float64, nr, 2) # allocate, need this?
                elk_apwfr[ia][l][io][:,:] = raw_apwfr[1:nr,1:2,io,l,ia]
            end
        end
    end

    return elk_apwfr
end

function compare_apwfr(apwfr, elk_apwfr)
    Natoms = size(apwfr, 1)
    lmaxapw = size(apwfr[1], 1) - 1
    Norder = size(apwfr[1][0], 1)
    for ia in 1:Natoms, l in 0:lmaxapw, io in 1:Norder
        mean_diff = mean(abs.(apwfr[ia][l][io] - elk_apwfr[ia][l][io]))
        println("ia=$ia l=$l io=$io mean_diff = $mean_diff")
    end
    return
end