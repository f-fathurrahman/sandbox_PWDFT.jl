using LinearAlgebra
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_forces_NN.jl")

function test_H2()
    atoms = Atoms(xyz_string="""
    3

    C   0.0   0.0    0.0
    O   1.5   0.0    0.0
    O   0.0   1.5    0.0
    """, in_bohr=true, LatVecs=gen_lattice_sc(16.0))
    atoms.Zvals = [4.0, 6.0]
    Zvals = atoms.Zvals
    Natoms = atoms.Natoms
    println(atoms)

    pw = PWGrid(15.0, atoms.LatVecs)

    F_NN = zeros(3,atoms.Natoms)
    calc_forces_NN!(pw, atoms, Zvals, F_NN)
    F_NN = 2*F_NN
    println("In Ry/au")
    println("Real space sum contribution: ")
    for ia in 1:Natoms
        @printf("%18.10f %18.10f %18.10f\n", F_NN[1,ia], F_NN[2,ia], F_NN[3,ia])
    end

    #F_NN = calc_forces_NN(atoms)
    ##display(F_NN'); println()
    #println("In Ry/au")
    #display(2*F_NN'); println()

end

test_H2()
