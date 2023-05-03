using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

function main()
    Ham, pwinput = init_Ham_from_pwinput()

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    nhm = Ham.pspotNL.nhm
    Natoms = Ham.atoms.Natoms
    Nspin = 1

    Nbecsum = Int64( nhm * (nhm + 1)/2 )
    println("Nbecsum = ", Nbecsum)

    becsum = zeros(Float64, Nbecsum, Natoms, Nspin)
    becsum[1:5,1,1] .= [1.0, 2.0, 3.0, 4.0, 5.0]
    becsum[1:5,2,1] .= [8.0, 8.0, 8.0, 8.0, 8.0]

    println("sum becsum before PAW_symmetrize: ", sum(becsum))
    for i in 1:5
        @printf("%3d %18.10f\n", i, becsum[i,1,1])
    end
    println()
    for i in 1:5
        @printf("%3d %18.10f\n", i, becsum[i,2,1])
    end

    PAW_symmetrize!(Ham.atoms, Ham.pspots, Ham.sym_info, Ham.pspotNL, becsum)

    println("sum becsum after PAW_symmetrize: ", sum(becsum))
    for i in 1:5
        @printf("%3d %18.10f\n", i, becsum[i,1,1])
    end
    println()
    for i in 1:5
        @printf("%3d %18.10f\n", i, becsum[i,2,1])
    end

end

main()