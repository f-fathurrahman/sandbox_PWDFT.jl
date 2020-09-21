using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("INC_gamma_only.jl")

function do_calc_forces(molname)

    Random.seed!(1234)

    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    
    ecutwfc = 15.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc )
    psis = randn_BlochWavefuncGamma(Ham)
    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )
    KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random,
        skip_initial_diag=true, NiterMax=200 )

    # Gamma-only
    F_Ps_loc  = zeros(3, Ham.atoms.Natoms)
    calc_forces_Ps_loc!(Ham, F_Ps_loc)
    println("F_Ps_loc:")
    display(F_Ps_loc'); println()

    F_Ps_loc_  = zeros(3,atoms.Natoms)
    calc_forces_Ps_loc!(Ham_, F_Ps_loc_)
    println("F_Ps_loc_:")
    display(F_Ps_loc_'); println()


    # Gamma-only
    F_Ps_nloc  = zeros(3, Ham.atoms.Natoms)
    calc_forces_Ps_nloc!(Ham, psis, F_Ps_nloc)
    println("F_Ps_nloc:")
    display(F_Ps_nloc'); println()

    F_Ps_nloc_  = zeros(3,atoms.Natoms)
    calc_forces_Ps_nloc!(Ham_, psiks, F_Ps_nloc_)
    println("F_Ps_nloc_:")
    display(F_Ps_nloc_'); println()

    forces = calc_forces(Ham, psis)
    forces_ = calc_forces(Ham_, psiks)

    println("forces = ")
    display(forces'); println()

    println("forces_ = ")
    display(forces_'); println()

end


function main()
    Nargs = length(ARGS)
    if Nargs >= 1
        molname = ARGS[1]
    else
        molname = "H2O"
    end
    do_calc_forces(molname)
end

main()
