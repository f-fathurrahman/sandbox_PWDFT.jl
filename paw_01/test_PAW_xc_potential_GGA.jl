using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

function main(;filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    becsum = PAW_atomic_becsum(Ham.atoms, Ham.pspots, Ham.pspotNL, Nspin=1)
    
    println("sum becsum before PAW_symmetrize: ", sum(becsum))
    PAW_symmetrize!(Ham, becsum)
    println("sum becsum after PAW_symmetrize: ", sum(becsum))

    atoms = Ham.atoms
    pspots = Ham.pspots
    pspotNL = Ham.pspotNL
    AE = true
    ia = 2
    isp = atoms.atm2species[ia]
    Nrmesh = pspots[isp].Nr
    l2 = (pspots[isp].lmax_rho + 1)^2
    Nspin = 1
    xc_calc = Ham.xc_calc

    # Calculate rho_lm
    rho_lm = zeros(Float64, Nrmesh, l2, Nspin)
    PAW_rho_lm!(AE, ia, atoms, pspots, pspotNL, becsum, rho_lm)
    println("sum rho_lm = ", sum(rho_lm))
    v_lm = zeros(Float64, Nrmesh, l2, Nspin)

    energy = PAW_xc_potential_GGA!( AE, ia,
        atoms, pspots, pspotNL, xc_calc,
        rho_lm, v_lm
    )

    println()
    println("energy for all ix = ", energy)
    println("sum abs vout_lm = ", sum(abs.(v_lm)))



end

main()