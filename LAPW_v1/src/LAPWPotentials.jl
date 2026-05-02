mutable struct LAPWPotentials
    vclmt::Vector{Vector{Float64}}
    vclir::Vector{Float64}
    epsxcmt::Vector{Vector{Float64}}
    epsxcir::Vector{Float64}
    vxcmt::Vector{Vector{Float64}}
    vxcir::Vector{Float64}
    bxcmt::Union{Nothing,Vector{Matrix{Float64}}} # ndmag
    bxcir::Union{Nothing,Matrix{Float64}} # ndmag
    #
    vsmt::Vector{Vector{Float64}}
    vsir::Vector{Float64}
    bsmt::Union{Nothing,Vector{Matrix{Float64}}} # ndmag
    bsir::Union{Nothing,Matrix{Float64}} # ndmag
    #
    vsig::Vector{ComplexF64}
end

function LAPWPotentials(
    atoms, mt_vars, pw, elk_input
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    npmt = mt_vars.npmt
    spinpol = elk_input.spinpol
    ndmag = elk_input.ndmag
    Npoints = prod(pw.Ns)

    # Hartree potential (MT and interstitial)
    vclmt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vclmt[ia] = zeros(Float64, npmt[isp])
    end
    vclir = zeros(Float64, Npoints)

    # Exchange correlation potentials
    epsxcmt = Vector{Vector{Float64}}(undef,Natoms)
    vxcmt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        epsxcmt[ia] = zeros(Float64, npmt[isp])
        vxcmt[ia] = zeros(Float64, npmt[isp])
    end
    if spinpol
        bxcmt = Vector{Matrix{Float64}}(undef,Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bxcmt[ia] = zeros(Float64, npmt[isp], ndmag)
        end
    else
        bxcmt = nothing
    end
    #
    epsxcir = zeros(Float64, Npoints)
    vxcir = zeros(Float64, Npoints)
    if spinpol
        bxcir = zeros(Float64, Npoints, ndmag)
    else
        bxcir = nothing
    end

    vsmt = Vector{Vector{Float64}}(undef, Natoms)
    vsir = zeros(Float64, Npoints)
    for ia in 1:Natoms
        isp = atm2species[ia]
        vsmt[ia] = zeros(Float64, npmt[isp])
    end

    if spinpol
        bsir = zeros(Float64, Npoints, ndmag)
        bsmt = Vector{Matrix{Float64}}(undef,Natoms)
        for ia in 1:Natoms
            isp = atm2species[ia]
            bsmt[ia] = zeros(Float64, npmt[isp], ndmag)
        end
    else
        bsir = nothing
        bsmt = nothing
    end

    vsig = zeros(ComplexF64, pw.gvec.Ng)

    return LAPWPotentials(
        vclmt, vclir,
        epsxcmt, epsxcir,
        vxcmt, vxcir,
        bxcmt, bxcir,
        vsmt, vsir,
        bsmt, bsir,
        vsig
    )
end