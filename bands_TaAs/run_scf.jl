using PWDFT
import Serialization

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main_scf()

    LatVecs = zeros(Float64, 3, 3)
    LatVecs[:,1] = [3.4696139977410336, 0.0000000000000000, 0.0000000000000000]
    LatVecs[:,2] = [0.0000000009593282, 3.4696139977410327, 0.0000000000000000]
    LatVecs[:,3] = [-1.7348069988705173, -1.7348069993501816, 5.8674449964211588]
    LatVecs *= ANG2BOHR

    atoms = Atoms(xyz_string=
        """
        4

        Ta       0.00000000       1.73480700       2.93297147
        Ta       0.00000000      -0.00000000       5.86669396
        As       1.73480700      -0.00000000       1.96869208
        As       1.73480700       1.73480700       4.90241458
        """,
        LatVecs = LatVecs
    )
    write_xsf("GEOM_exported.xsf", atoms)

    pspfiles = [
        joinpath(DIR_PSP, "Ta-q5.gth"),
        joinpath(DIR_PSP, "As-q5.gth")
    ]
    ecutwfc = 20.0

    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc, meshk=[6,6,3], extra_states=5
    )
    KS_solve_SCF!( Ham, mix_method = "rpulay", use_smearing = true)

    Serialization.serialize("Ham.data", Ham)

end

#main_scf()
