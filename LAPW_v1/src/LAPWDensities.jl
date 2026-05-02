mutable struct LAPWDensities
    rhomt::Vector{Vector{Float64}}
    rhoir::Vector{Float64}
    magmt::Union{Nothing,Vector{Matrix{Float64}}}
    magir::Union{Nothing,Matrix{Float64}}
end

# Initialize densities by callin rhoinit and maginit
function LAPWDensities(
    atoms, atsp_vars, mt_vars, pw, elk_input
)
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    npmt = mt_vars.npmt
    Npoints = prod(pw.Ns)
    spinpol = elk_input.spinpol

    # Initialize rhomt and rhoir
    rhomt = Vector{Vector{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        rhomt[ia] = zeros(Float64, npmt[isp])
    end
    rhoir = zeros(Float64, Npoints)
    #
    rhoinit!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir )

    # magnetization
    if spinpol
        ncmag = elk_input.ncmag
        ndmag = elk_input.ndmag
        bfcmt = elk_input.bfcmt0
        bfieldc = elk_input.bfieldc
        magmt = Vector{Matrix{Float64}}(undef, Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            magmt[ia] = zeros(Float64, npmt[isp], ndmag)
        end
        magir = zeros(Float64, Npoints, ndmag)
        maginit!(atoms, rhomt, rhoir, ncmag, ndmag, bfcmt, bfieldc, magmt, magir)
    else
        magmt = nothing
        magir = nothing
    end

    return LAPWDensities(rhomt, rhoir, magmt, magir)
end