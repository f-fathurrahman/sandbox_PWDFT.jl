using Printf, LinearAlgebra, Serialization
using FFTW
using SpecialFunctions: erf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("calc_stress_hartree.jl")
include("eval_dVloc_G.jl")
include("calc_stress_Ps_loc.jl")
include("calc_stress_xc.jl")
include("calc_stress_nlcc.jl")
include("calc_stress_kinetic.jl")
include("calc_Deff.jl")
include("gen_us_dj.jl")
include("gen_us_dy.jl")
include("calc_stress_Ps_nloc.jl")

function main()
    Ham = Serialization.deserialize("Hamiltonian.dat")
    psiks = Serialization.deserialize("psiks.dat")

#=
    println(Ham.energies)

    stress_hartree = zeros(Float64, 3, 3)
    calc_stress_hartree!( Ham.pw, Ham.rhoe, Ham.energies.Hartree, stress_hartree )
    stress_hartree *= 2.0 # convert to Ry
    println("\nstress_hartree (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_hartree[i,1], stress_hartree[i,2], stress_hartree[i,3])
    end


    stress_Ps_loc = zeros(Float64, 3, 3)
    calc_stress_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, stress_Ps_loc )
    stress_Ps_loc *= 2.0 # convert to Ry
    println("\nstress_Ps_loc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_Ps_loc[i,1], stress_Ps_loc[i,2], stress_Ps_loc[i,3])
    end


    stress_xc = zeros(Float64, 3, 3)
    calc_stress_xc!( Ham.pw, Ham.potentials, Ham.rhoe, Ham.energies.XC, stress_xc )
    stress_xc *= 2.0 # convert to Ry
    println("\nstress_xc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_xc[i,1], stress_xc[i,2], stress_xc[i,3])
    end


    stress_nlcc = zeros(Float64, 3, 3)
    calc_stress_nlcc!( Ham.atoms, Ham.pspots, Ham.pw, Ham.xc_calc, Ham.xcfunc,
        Ham.rhoe, Ham.rhoe_core, stress_nlcc )
    stress_nlcc *= 2.0 # convert to Ry
    println("\nstress_nlcc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_nlcc[i,1], stress_nlcc[i,2], stress_nlcc[i,3])
    end


    stress_kin = zeros(Float64, 3, 3)
    calc_stress_kinetic!( Ham.pw, Ham.electrons, psiks, stress_kin )
    stress_kin *= 2.0 # convert to Ry
    println("\nstress_kin (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_kin[i,1], stress_kin[i,2], stress_kin[i,3])
    end
=#

    stress_Ps_nloc = zeros(Float64, 3, 3)
    calc_stress_Ps_nloc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.pspotNL, Ham.electrons,
        psiks, stress_Ps_nloc )
    stress_Ps_nloc *= 2.0 # convert to Ry
    println("\nstress_Ps_nloc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", stress_Ps_nloc[i,1], stress_Ps_nloc[i,2], stress_Ps_nloc[i,3])
    end

    return
end

main()
