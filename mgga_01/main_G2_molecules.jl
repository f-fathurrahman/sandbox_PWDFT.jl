using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_PSP_PBE = joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("../get_default_psp.jl")

function do_calc(molname)

    Random.seed!(1234)

    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms, xcfunc="PBE")
    
    ecutwfc = 15.0
    
    Ham = Hamiltonian(atoms, pspfiles, ecutwfc, xcfunc="SCAN", use_symmetry=false)
    println(Ham.pw)
    psiks = rand_BlochWavefunc(Ham)

    KS_solve_Emin_PCG!(Ham, psiks, print_final_ebands=true)
    #KS_solve_SCF!(Ham, psiks, mix_method="rpulay")
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
