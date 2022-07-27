using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main()

    Random.seed!(1234)

    # Atoms
    atoms = Atoms(xyz_string="""
    2

    Si  3.505  2.5   2.5
    Si  0.0    2.5   2.5
    """, in_bohr=true, LatVecs = gen_lattice_sc(10.0) )

    # Initialize Hamiltonian
    pspfiles = ["Si-q4_mod.gth"]
    #pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 10.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false, extra_states=4 )
    println(Ham)

    KS_solve_SCF!( Ham, mix_method="pulay", use_smearing=true )
end

main()
