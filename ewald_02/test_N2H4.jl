using LinearAlgebra
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_forces_NN.jl")

function main()
    # Coords are in angstrom
    atoms = Atoms(xyz_string="""
    6
    
    N       5.94821400       6.81171100       5.22639100
    N       5.94821400       5.37379300       5.22639100
    H       6.15929600       7.18550400       6.15196500
    H       5.00000000       7.09777800       5.00000000
    H       5.73713200       5.00000000       6.15196500
    H       6.89642800       5.08772600       5.00000000
    """, LatVecs=gen_lattice_sc(16.0))
    atoms.Zvals = [5.0, 1.0]
    Zvals = atoms.Zvals
    Natoms = atoms.Natoms
    println(atoms)

    pw = PWGrid(15.0, atoms.LatVecs)

    F_NN = zeros(3,atoms.Natoms)
    calc_forces_NN!(pw, atoms, Zvals, F_NN)
    fill!(F_NN, 0.0)
    @time calc_forces_NN!(pw, atoms, Zvals, F_NN)
    F_NN = 2*F_NN
    println("In Ry/au")
    println("Real space sum contribution: ")
    for ia in 1:Natoms
        @printf("%18.10f %18.10f %18.10f\n", F_NN[1,ia], F_NN[2,ia], F_NN[3,ia])
    end

    F_NN = calc_forces_NN(atoms)
    @time F_NN = calc_forces_NN(atoms)
    #display(F_NN'); println()
    println("In Ry/au")
    display(2*F_NN'); println()

end

main()
