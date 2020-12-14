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

function create_Ham_H2()
    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP_PBE, "H-q1.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="SCAN", use_symmetry=false )
end

function main()

    Random.seed!(1234)
    
    Ham = create_Ham_H2()
    println(Ham.pw)
    psiks = rand_BlochWavefunc(Ham)

    KS_solve_Emin_PCG!(Ham, psiks)
    #KS_solve_SCF!(Ham, psiks, mix_method="rpulay", print_final_ebands=true)
end

main()
