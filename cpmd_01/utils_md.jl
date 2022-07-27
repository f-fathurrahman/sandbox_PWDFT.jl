function write_traj_files(
    filetraj, fileetot,
    atoms, forces,
    Etot_conserved, Etot, Ekin_ions, Ekin_elec,
    dt, iter
)
    # Forces are in Ha/bohr
    # Xcrysden expect the forces are in Ha/angstrom
    F = forces*(1.0/BOHR2ANG)
    #
    @printf(filetraj, "%d  Etot_conserved = %18.10f\n\n", atoms.Natoms, Etot_conserved)
    @printf(fileetot, "%18.10f %18.10f %18.10f %18.10f %18.10f\n",
            AU_PS*iter*dt, Etot_conserved, Etot, Ekin_ions, Ekin_elec)
    for ia in 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        r = atoms.positions*BOHR2ANG # convert from bohr to angstrom
        @printf(filetraj, "%3s %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                atoms.SpeciesSymbols[isp],
                r[1,ia], r[2,ia], r[3,ia],
                F[1,ia], F[2,ia], F[3,ia]
        ) # printf
    end
    flush(filetraj)
    flush(fileetot)

    return
end


function info_md(atoms, energies, forces, iter)
    @printf("\nMD Iter = %3d, Etot = %18.10f\n", iter, sum(energies))
    println("Forces = ")
    for ia in 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        atsymb = atoms.SpeciesSymbols[isp]
        @printf("%3s %18.10f %18.10f %18.10f\n", atsymb,
            forces[1,ia], forces[2,ia], forces[3,ia])
    end
end