using Printf, LinearAlgebra, Serialization
using FFTW
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("calc_stress_hartree.jl")

function main()
    Ham = Serialization.deserialize("Hamiltonian.dat")
    psiks = Serialization.deserialize("psiks.dat")

    println(Ham.energies)

    stress_hartree = zeros(Float64, 3, 3)
    calc_stress_hartree!( Ham.pw, Ham.rhoe, Ham.energies.Hartree, stress_hartree )

    stress_hartree *= 2.0 # convert to Ry
    println("stress_hartree (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_hartree[i,1], stress_hartree[i,2], stress_hartree[i,3])
    end
    display(stress_hartree); println()

    return
end

main()
