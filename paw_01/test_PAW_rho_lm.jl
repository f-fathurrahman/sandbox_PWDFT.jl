using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("PAW_atomic_becsum.jl")
include("PAW_symmetrize.jl")
include("PAW_rho_lm.jl")

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
    isp = Ham.atoms.atm2species[ia]
    Nrmesh = Ham.pspots[isp].Nr
    l2 = (Ham.pspots[isp].lmax_rho + 1)^2

    rho_lm = zeros(Float64, Nrmesh, l2, Nspin)
    AE = true
    PAW_rho_lm!(AE, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm)
    println("sum rho_lm = ", sum(rho_lm))

end

main()