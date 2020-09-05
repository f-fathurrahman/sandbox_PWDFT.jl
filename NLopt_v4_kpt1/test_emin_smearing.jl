using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")
include("smearing.jl")
include("SubspaceRotations.jl")
include("ElecGradient.jl")
include("ElecVars.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")
include("setup_guess_wavefunc.jl")
include("KS_solve_Emin_SD_Haux.jl")

Random.seed!(1234)

function main(Ham::Hamiltonian; kT=0.01)
    psiks = rand_BlochWavefunc(Ham)
    #setup_guess_wavefunc!( Ham, psiks, startingrhoe=:gaussian, skip_initial_diag=false )
    setup_guess_wavefunc!( Ham, psiks, startingrhoe=:random, skip_initial_diag=false )
    evars = ElecVars(Ham, psiks)    

    Ham.electrons.ebands[:] = evars.Hsub_eigs[:] # Initialize Haux to eigvals of Hsub
    
    println("Initial guess:")
    print_ebands_Hsub_eigs(Ham, evars)

    KS_solve_Emin_SD_Haux!( Ham, evars, NiterMax=200, kT=kT )
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)
end

function main_scf(Ham::Hamiltonian; kT=0.01)
    #KS_solve_SCF!( Ham, mix_method="anderson", use_smearing=true, kT=kT)
    KS_solve_SCF_potmix!( Ham, use_smearing=true, kT=kT)
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)
end

#main(create_Ham_O2(), kT=0.01)
main_scf(create_Ham_O2(), kT=0.01)

#main(create_Ham_atom("Ni", "Ni-q10.gth"), kT=0.01)
#main_scf(create_Ham_atom("Ni", "Ni-q10.gth"), kT=0.01)

#main(create_Ham_atom("Co", "Co-q9.gth"), kT=0.01)
#main_scf(create_Ham_atom("Co", "Co-q9.gth"), kT=0.01)

#main(create_Ham_atom_Al_smearing(), kT=0.001)
#main_scf(create_Ham_atom_Al_smearing(), kT=0.001)

#main(create_Ham_atom_C_smearing(), kT=0.001)
#main_scf(create_Ham_atom_C_smearing(), kT=0.001)

#main(create_Ham_atom_Si_smearing(), kT=0.001)
#main_scf(create_Ham_atom_Si_smearing(), kT=0.001)

#main(create_Ham_atom_Pt_smearing(), kT=0.001)
#main_scf(create_Ham_atom_Pt_smearing(), kT=0.001)