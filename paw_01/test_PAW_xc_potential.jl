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

    Nspin = 1
    ia = 1
    AE = true

    isp = Ham.atoms.atm2species[ia]
    Nrmesh = Ham.pspots[isp].Nr
    l2 = (Ham.pspots[isp].lmax_rho + 1)^2

    rho_lm = zeros(Float64, Nrmesh, l2, Nspin)
    PAW_rho_lm!(AE, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm)
    println("sum rho_lm = ", sum(rho_lm))

    v_lm = zeros(Float64, Nrmesh, l2, Nspin)
    @views energy = PAW_xc_potential!(
        AE, ia,
        Ham.atoms, Ham.pspots, Ham.pspotNL, Ham.xc_calc,
        rho_lm, v_lm
    )

    println("\nEnd result:")
    println("After PAW_xc_potential: ")
    println("sum v_lm = ", sum(v_lm))
    println("energy = ", energy)

end

main()