using Printf
using LinearAlgebra
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("gen_neighbor_shells.jl")
include("calc_stress_NN.jl")

function main()
    Ham, pwinput = init_Ham_from_pwinput()

    stress_NN = zeros(Float64, 3, 3)
    calc_stress_NN!( Ham.atoms, Ham.pw, stress_NN )

    stress_NN *= 2.0 # convert to Ry
    println("stress_NN (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_NN[i,1], stress_NN[i,2], stress_NN[i,3])
    end
    display(stress_NN); println()

    return
end

main()
