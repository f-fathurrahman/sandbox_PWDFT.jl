push!(LOAD_PATH, pwd())

using PWDFT
using LAPWDFT

function create_atomic_vars(lattice_vars)
    maxatoms = 200
    maxspecies = 8
    atposl = zeros(3,maxatoms,maxspecies)
    Nspecies = 2
    Natoms = [2,1]
    # species 1, atom 1 and 2
    atposl[:,1,1] = [0.0, 0.0, 0.0]
    atposl[:,2,1] = [0.25, 0.25, 0.25]
    # species 2, atom 1
    atposl[:,1,2] = [0.35, 0.35, 0.15]
    atomic_vars = AtomicVars(Nspecies, Natoms, atposl, lattice_vars)
end

function main()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        3

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        Pt  0.35  0.35  0.15
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    println(atoms)

    atomic_vars = create_atomic_vars(LatticeVars(atoms.LatVecs))

    Natoms = atomic_vars.Natoms
    Nspecies = atomic_vars.Nspecies
    for isp in 1:Nspecies
        println("atposl")
        display(atomic_vars.atposl[:,1:Natoms[isp],isp]); println()
        println("atposc")
        display(atomic_vars.atposc[:,1:Natoms[isp],isp]); println()
    end
end

main()