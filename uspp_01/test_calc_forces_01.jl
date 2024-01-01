using Printf
import Serialization
using FFTW

using PWDFT

include("calc_forces_scf_corr.jl")
include("calc_forces_Ps_nloc.jl")

function main()
    
    Ham = Serialization.deserialize("Hamiltonian.dat")
    psiks = Serialization.deserialize("psiks.dat")

    atoms = Ham.atoms
    pw = Ham.pw
    pspots = Ham.pspots
    electrons = Ham.electrons
    pspotNL = Ham.pspotNL
    atsymbs = atoms.atsymbs
    Natoms = atoms.Natoms
    potentials = Ham.potentials

    #
    # Nonlocal pspot components
    #
    F_Ps_nloc = zeros(Float64,3,atoms.Natoms)
    my_calc_forces_Ps_nloc!(
        atoms, pw, pspots,
        electrons, pspotNL, potentials, psiks, F_Ps_nloc
    )

    F_Ps_nloc[:] .= F_Ps_nloc[:]*2.0
    println("F_Ps_nloc: (in Ry/bohr)")
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end

    #
    # SCF correction
    #
    F_scf_corr = zeros(Float64,3,atoms.Natoms)
    my_calc_forces_scf_corr!(
        atoms, pw, pspots, potentials, F_scf_corr
    )
    F_scf_corr[:] .*= 2.0 # convert to Ry/bohr
    println("F_scf_corr: (in Ry/bohr)")
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_scf_corr[1,ia], F_scf_corr[2,ia], F_scf_corr[3,ia] )
    end

end

main()
