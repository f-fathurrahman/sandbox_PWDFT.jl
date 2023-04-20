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
include("../ld1/RadialGrid.jl")
include("radial_hartree.jl")
include("PAW_rho_lm.jl")
include("PAW_lm2rad.jl")
include("PAW_rad2lm.jl")
include("PAW_h_potential.jl")
include("PAW_xc_potential.jl")
include("PAW_potential.jl")


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
    println("EHxc_paw = ", EHxc_paw)

end

main()