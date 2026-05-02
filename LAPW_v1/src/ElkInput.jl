struct ElkInput
    LatVecs::Matrix{Float64}
    Nspecies::Int64
    species_files::Vector{String}
    natoms_per_species::Vector{Int64}
    atomic_positions::Vector{Matrix{Float64}}
    is_molecule::Bool
    ngridk::Vector{Int64}
    nempty::Union{Int64,Nothing}
    spinpol::Bool
    bfieldc::Vector{Float64}
    spinorb::Bool
    cmagz::Bool
    nosource::Bool
    spinsprl::Bool
    lradstp::Int64
    ncmag::Bool
    rgkmax::Float64
    gmaxvr::Float64
    bfcmt0::Matrix{Float64}
    epslat::Float64
    ndmag::Int64
end
# NOTE: species_files are assumed to be located in the current directory
# For bfieldc, I think it is more flexible to use Vector{Float64} rather than Tuple
#
# bfcmt0 is already "flattened"

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
    nempty = nothing
    spinpol = false
    scale = 1.0
    bfieldc = zeros(Float64, 3)
    spinorb = false
    cmagz = false
    nosource = false
    spinsprl = false
    lradstp = 4
    ncmag = false
    rgkmax = 7.0
    gmaxvr = 12.0
    epslat = 1e-6
    ndmag = -1 # some invalid value

    # Need to total number of
    #bfcmt = zeros(Float64, 3, atoms.Natoms)
    #bfcmt0 = zeros(Float64, 3, atoms.Natoms)
    bfcmt0_in = Vector{Matrix{Float64}}()

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
                bfcmt0_sp = zeros(Float64, 3, natmsp)
                # Start reading atomic positions
                for ia in 1:natmsp
                    iline += 1
                    ll = split(lines[iline], " ", keepempty=false)
                    @info "ll = $(ll)"
                    for i in 1:3
                        atpos_sp[i,ia] = parse(Float64, ll[i])
                    end
                    # This is bfcmt0
                    if length(ll) > 3
                        for i in 1:3
                            bfcmt0_sp[i,ia] = parse(Float64, ll[3+i])
                        end
                    end
                end
                push!(atomic_positions, atpos_sp)
                push!(bfcmt0_in, bfcmt0_sp)
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
        if l == "nempty"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            nempty = parse(Int64, ll)
        end
        #
        if l == "lradstp"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            lradstp = parse(Int64, ll)
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
        if l == "ncmag"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                ncmag = true
            end 
        end
        #
        if l == "rgkmax"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            rgkmax = parse(Float64, ll)
        end
        #
        if l == "gmaxvr"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            gmaxvr = parse(Float64, ll)
        end
        #
        if l == "scale"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            scale = parse(Float64, ll)
        end
        #
        if l == "bfieldc"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)
            bfieldc[1] = parse(Float64, ll[1])
            bfieldc[2] = parse(Float64, ll[2])
            bfieldc[3] = parse(Float64, ll[3])
        end
        #
        if l == "spinorb"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                spinorb = true
            end
        end
        #
        if l == "cmagz"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                cmagz = true
            end
        end
        #
        if l == "nosource"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                nosource = true
            end
        end
        #
        if l == "spinsprl"
            iline += 1
            ll = split(lines[iline], " ", keepempty=false)[1]
            if lowercase(ll) == ".true."
                spinsprl = true
            end
        end

    end

    # Scale LatVecs
    LatVecs *= scale

    # "Flatten" bfcmt0
    Natoms = sum(natoms_per_species)
    bfcmt0 = zeros(Float64, 3, Natoms)
    ia = 0
    for isp in 1:Nspecies
        for ias in 1:natoms_per_species[isp]
            ia += 1
            bfcmt0[:,ia] = bfcmt0_in[isp][:,ias]
        end
    end

    if is_molecule
        ngridk = ones(Int64, 3)
        # This is the default anyway
    end

    bfieldc0 = bfieldc # need these two?
    # check for collinearity in the z-direction and set the dimension of the
    # magnetization and exchange-correlation vector fields
    if spinpol
        ndmag = 1
        if ( abs(bfieldc0[1]) > epslat) || (abs(bfieldc0[2]) > epslat)
            ndmag = 3
            @info "ndmag is set to 3"
        end
        for ia in 1:Natoms
            if (abs(bfcmt0[1,ia]) > epslat || (abs(bfcmt0[2,ia]) > epslat))
                ndmag = 3
                @info "ndmag is set to 3"
            end
        end
        # spin-orbit coupling is non-collinear in general
        if spinorb
            ndmag = 3
        end
        # source-free fields and spin-spirals must be non-collinear
        if nosource || spinsprl
            ndmag = 3
            cmagz = false
        end
        # force collinear magnetism along the z-axis if required
        if cmagz
            ndmag = 1
        end
    else
        ndmag = 0
    end

    @assert ndmag >= 0

    return ElkInput(
        LatVecs, Nspecies, species_files,
        natoms_per_species, atomic_positions,
        is_molecule,
        ngridk, nempty,
        spinpol, bfieldc,
        spinorb, cmagz, nosource, spinsprl,
        lradstp, ncmag, rgkmax, gmaxvr,
        bfcmt0, epslat, ndmag
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
