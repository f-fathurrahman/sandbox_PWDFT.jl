#
# A naive function to read pwscf input
# Only a subset of possible combinations of parameters are considered.
#
function read_pwscf_input( filename::String )

    # Default values, some are not valid for PWSCF
    acell = -1.0

    Natoms = 0
    Nspecies = 0

    LatVecs = zeros(3,3)
    is_parse_cell = false
    N_parse_cell = 0

    xyz_string_frac = ""
    xyz_string_angstrom = ""
    xyz_string_bohr = ""
    is_parse_xyz = false
    N_parse_xyz = 0

    is_parse_kpoints = false
    N_parse_kpoints = 0
    meshk = [0, 0, 0]

    in_angstrom = false
    in_fraction = false
    in_bohr = false

    f = open(filename, "r")
    
    while !eof(f)
        
        l = readline(f)
        
        if occursin("  A =", l)
            ll = split(l, "=")
            acell = parse(Float64,ll[end])*ANG2BOHR
        end

        # Read number of atoms
        if occursin("nat =", l)
            ll = split(l, "=", keepempty=false)
            Natoms = parse(Int64, ll[end])
            println("Read Natoms = ", Natoms)
        end

        # Read number of species (ntyp)
        if occursin("ntyp =", l)
            ll = split(l, "=", keepempty=false)
            Nspecies = parse(Int64, ll[end])
            println("Read Nspecies = ", Nspecies)
        end

        # Read cell parameters
        # FIXME: They are assumed to be given in bohr !!!!
        if occursin("CELL_PARAMETERS", l)
            is_parse_cell = true
        end

        if is_parse_cell && N_parse_cell <= 3
            if N_parse_cell == 0
                N_parse_cell = N_parse_cell + 1
                continue
            end
            ll = split(l, " ", keepempty=false)
            LatVecs[1,N_parse_cell] = parse(Float64, ll[1])
            LatVecs[2,N_parse_cell] = parse(Float64, ll[2])
            LatVecs[3,N_parse_cell] = parse(Float64, ll[3])
            N_parse_cell = N_parse_cell + 1
        end


        # XXX: Natoms should be read before
        if occursin("ATOMIC_POSITIONS", l)
            is_parse_xyz = true
            if occursin("crystal", l)
                in_fraction = true
                xyz_string_frac = xyz_string_frac*string(Natoms)*"\n\n"
            end
            if occursin("angstrom", l)
                in_angstrom = true
                xyz_string_angstrom = xyz_string_angstrom*string(Natoms)*"\n\n"
            end
            if occursin("bohr", l)
                in_bohr = true
                xyz_string_bohr = xyz_string_bohr*string(Natoms)*"\n\n"
            end
        end

        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_fraction
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_frac = xyz_string_frac*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end

        # Atomic positions are given in angstrom
        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_angstrom
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_angstrom = xyz_string_angstrom*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end

        # Atomic positions are given in bohr
        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_bohr
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_bohr = xyz_string_bohr*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end


        if occursin("K_POINTS", l)
            is_parse_kpoints = true
        end

        if is_parse_kpoints && N_parse_kpoints <= 1
            if N_parse_kpoints == 0
                N_parse_kpoints = N_parse_kpoints + 1
                continue
            end
            ll = split(l, " ", keepempty=false)
            meshk[1] = parse(Int64,ll[1])
            meshk[2] = parse(Int64,ll[2])
            meshk[3] = parse(Int64,ll[3])
            N_parse_kpoints = N_parse_kpoints + 1
        end

    end
    close(f)

    # Only scale LatVecs with acell if it is a positive value
    if acell > 0.0
        LatVecs = acell*LatVecs
    end

    if in_fraction
        println(xyz_string_frac)
        atoms = init_atoms_xyz_string( xyz_string_frac, in_bohr=true )
        atoms.positions = LatVecs*atoms.positions # convert to bohr
    #
    elseif in_angstrom
        println(xyz_string_angstrom)
        atoms = init_atoms_xyz_string( xyz_string_angstrom )
    #
    elseif in_bohr
        println(xyz_string_bohr)
        atoms = init_atoms_xyz_string( xyz_string_bohr, in_bohr=true )
    else
        error("Cannot read atomic positions")
    end
    # Set unit lattice vectors manually
    atoms.LatVecs = LatVecs

    return atoms, meshk
end