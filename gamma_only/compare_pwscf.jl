using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

include("INC_gamma_only.jl")

function run_pwscf(Ham; etot_conv_thr=2e-6)
    run(`rm -rfv TEMP_pwscf/\*`)
    write_pwscf(Ham, prefix_dir="TEMP_pwscf",
        gamma_only=true, etot_conv_thr=etot_conv_thr,
        nosym=true, mixing_beta=0.1)
    cd("./TEMP_pwscf")
    run(pipeline(`pw.x`, stdin="PWINPUT", stdout="LOG1"))
    cd("../")

    pwscf_energies = read_pwscf_etotal("TEMP_pwscf/LOG1")
    pwscf_forces = read_pwscf_forces("TEMP_pwscf/LOG1")
    return pwscf_energies, pwscf_forces
end



function do_calc(molname)

    Random.seed!(1234)

    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    
    ecutwfc = 15.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc)
    psis = randn_BlochWavefuncGamma(Ham)
    
    KS_solve_Emin_PCG_dot!(Ham, psis, NiterMax=200, etot_conv_thr=1e-8)
    energies = Ham.energies
    forces = calc_forces(Ham, psis)

    pwscf_energies, pwscf_forces = run_pwscf(Ham, etot_conv_thr=2e-8)

    println("PWSCF energies: ")
    println(pwscf_energies)

    display(forces'); println()
    display(pwscf_forces'); println()
end


function main()
    Nargs = length(ARGS)
    if Nargs >= 1
        molname = ARGS[1]
    else
        molname = "H2O"
    end
    do_calc(molname)
end

main()
