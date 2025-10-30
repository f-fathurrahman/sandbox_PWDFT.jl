struct ElkInput
    LatVecs::Matrix{Float64}
    Nspecies::Int64
    species_files::Vector{String}
    natoms_per_species::Vector{Int64}
    atomic_positions::Vector{Matrix{Float64}}
    is_molecule::Bool
    ngridk::Vector{Int64}
    spinpol::Bool
end
# NOTE: species_files are assumed to be located in the current directory


function read_elk_input()
    f = open("elk.in", "r")
    lines = readlines(f)
    close(f)

    # Some default values
    LatVecs = zeros(Float64, 3, 3)
    Nspecies = 0
    species_files = Vector{String}()
    natoms_per_species = Vector{Int64}()
    atomic_positions = Vector{Matrix{Float64}}()
    is_molecule = false
    ngridk = ones(Int64, 3)
    spinpol = false
    scale = 1.0

    Nlines = length(lines)
    iline = 0
    while true
        #
        iline += 1
        if iline >= Nlines
            break
        end
        #
        l = lines[iline]
        #
        if l == "avec"
            @info "Processing avec"
            iline += 1; l1 = lines[iline]
            iline += 1; l2 = lines[iline]
            iline += 1; l3 = lines[iline]
            # Read LatVecs by column
            # The 1st vector
            ll = split(l1, " ", keepempty=false)
            for i in 1:3
                LatVecs[i,1] = parse(Float64, ll[i])
            end
            # The 2nd vector
            ll = split(l2, " ", keepempty=false)
            for i in 1:3
                LatVecs[i,2] = parse(Float64, ll[i])
            end
            # The 3rd vector
            ll = split(l3, " ", keepempty=false)
            for i in 1:3
                LatVecs[i,3] = parse(Float64, ll[i])
            end
            @info "End of processing avec"
        end
        #
        if l == "atoms"
            @info "Processing atoms"
            # Read number of species
            iline += 1
            l = split(lines[iline], " ", keepempty=false)[1]
            Nspecies = parse(Int64, l)
            @info "Nspecies = $(Nspecies)"
            #
            for isp in 1:Nspecies
                iline += 1
                sp_path = split(lines[iline], " ", keepempty=false)[1]
                # Remove '
                sp_path = replace(sp_path, "'" => "")
                push!(species_files, sp_path)
                #
                iline += 1
                ll = split(lines[iline], " ", keepempty=false)[1]
                natmsp = parse(Int64, ll)
                @info "natmsp = $(natmsp)"
                push!(natoms_per_species, natmsp)
                atpos_sp = zeros(Float64, 3, natmsp)
                # Start reading atomic positions
                for ia in 1:natmsp
                    iline += 1
                    ll = split(lines[iline], " ", keepempty=false)
                    @info "ll = $(ll)"
                    for i in 1:3
                        atpos_sp[i,ia] = parse(Float64, ll[i])
                    end
                end
                push!(atomic_positions, atpos_sp)
            end 
            @info "End of processing atoms"
        end
        #
        if l == "molecule"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                is_molecule = true
            end
        end
        #
        if l == "ngridk"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)
            for i in 1:3
                ngridk[i] = parse(Int64, ll[i])
            end
        end
        #
        if l == "spinpol"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                spinpol = true
            end 
        end
        #
        if l == "scale"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            scale = parse(Float64, ll)
        end
    end

    # Scale LatVecs
    LatVecs *= scale

    return ElkInput(
        LatVecs, Nspecies, species_files,
        natoms_per_species, atomic_positions,
        is_molecule,
        ngridk, spinpol
    )
end


function create_atoms_from_elk_input(elk_input)
    Natoms = sum(elk_input.natoms_per_species)
    Nspecies = elk_input.Nspecies
    # Build atoms string
    atoms_string = "$Natoms\n\n"
    for isp in 1:Nspecies
        species_name = replace(elk_input.species_files[isp], ".in" => "")
        if length(species_name) > 3
            @warn "Too loop species_name: $(species_name)"
        end
        atpos = elk_input.atomic_positions[isp]
        for ias in 1:elk_input.natoms_per_species[isp]
            atoms_string *= "$species_name $(atpos[1,ias]) $(atpos[2,ias]) $(atpos[3,ias])\n"
        end
    end
    print(atoms_string)
    if elk_input.is_molecule
        return Atoms(xyz_string=atoms_string, in_bohr=true, LatVecs=elk_input.LatVecs)
    else
        # The coordinates are fractional
        return Atoms(xyz_string_frac=atoms_string, in_bohr=true, LatVecs=elk_input.LatVecs)
    end
end
