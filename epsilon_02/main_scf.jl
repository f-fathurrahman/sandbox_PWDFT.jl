using Printf
using PWDFT
import Serialization

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function init_Hamiltonian()    
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 20.0
    return Hamiltonian( atoms, pspfiles, ecutwfc,
        xcfunc="VWN", meshk=[6,6,6], shiftk=[1,1,1] )
end


function main_scf()
    Random.seed!(1234)
    Ham = init_Hamiltonian()
    println(Ham)
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.5 )
    Serialization.serialize("Ham.data", Ham)
end

main_scf()
