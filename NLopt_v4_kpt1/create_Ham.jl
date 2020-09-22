function create_Ham_O2()
    atoms = Atoms(xyz_string=
        """
        2
        
        O       0.622978       0.00000000       0.00000000
        O      -0.622978       0.00000000       0.00000000
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

function create_Ham_Si_dimer()
    atoms = Atoms(xyz_string=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.0   0.0   5.5
        """, in_bohr=true, LatVecs=diagm([10.0, 10.0, 15.0])
    )
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end


function create_Ham_atom(atsymb::String, pspfile::String)
    atoms = Atoms(xyz_string_frac=
        """
        1

        $atsymb  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, pspfile)]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

function create_Ham_atom_Pt_smearing(; a=16.0)
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(a))
    pspfiles = [joinpath(DIR_PSP, "Pt-q10.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

# without extra_states
function create_Ham_atom_Pt()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "Pt-q10.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function create_Ham_atom_Al_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

function create_Ham_atom_C_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        C  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "C-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end

function create_Ham_atom_Si_smearing()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Si  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, extra_states=4 )
end