using Printf
using OffsetArrays
using LinearAlgebra
using Serialization

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

    ddd_paw = Ham.pspotNL.paw.ddd_paw
    e_cmp = zeros(Float64, Ham.atoms.Natoms, 2, 2)
    EHxc_paw = PAW_potential!( Ham.atoms, Ham.pspots, Ham.pspotNL, Ham.xc_calc,
        becsum, ddd_paw, e_cmp
    )

    serialize("ddd_paw_jl.dat", ddd_paw)
    println("EHxc_paw = ", EHxc_paw)
    println("sum ddd_paw = ", sum(ddd_paw))

end

main(filename="PWINPUT")