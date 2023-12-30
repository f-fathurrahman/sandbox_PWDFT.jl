using Printf
import Serialization

using LinearAlgebra
using FFTW
using PWDFT

include("calc_forces_scf_corr.jl")

function main()
    
    Ham = Serialization.deserialize("Hamiltonian.dat")
    psiks = Serialization.deserialize("psiks.dat")

    atoms = Ham.atoms
    pw = Ham.pw
    pspots = Ham.pspots
    potentials = Ham.potentials
    F_scf_corr = zeros(Float64,3,atoms.Natoms)
    my_calc_forces_scf_corr!(
        atoms, pw, pspots, potentials, F_scf_corr
    )
    atsymbs = atoms.atsymbs
    Natoms = atoms.Natoms
    F_scf_corr[:] .*= 2.0 # convert to Ry/bohr
    println("F_scf_corr: (in Ry/bohr)")
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_scf_corr[1,ia], F_scf_corr[2,ia], F_scf_corr[3,ia] )
    end

end

main()
