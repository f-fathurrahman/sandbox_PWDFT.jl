struct POSCAR
    lattice::Matrix{Float64}
    symbols::Vector{String}
    coords::Matrix{Float64}
    coordinate_type::String
end


function read_poscar(filename::String)
    lines = strip.(readlines(filename))

    idx = 1

    comment = lines[idx]
    idx += 1

    scale = parse(Float64, lines[idx])
    idx += 1

    lattice = zeros(3,3)

    for i in 1:3
        lattice[i,:] .= parse.(Float64, split(lines[idx]))
        idx += 1
    end

    lattice .*= scale

    tokens = split(lines[idx])

    # Determine whether this line contains symbols or counts
    has_symbols = any(x -> occursin(r"[A-Za-z]", x), tokens)

    if has_symbols
        species = tokens
        idx += 1

        counts = parse.(Int, split(lines[idx]))
        idx += 1
    else
        error("POSCAR without species names is not supported")
    end

    # Expand symbols
    symbols = String[]
    for (sp,n) in zip(species, counts)
        append!(symbols, fill(sp,n))
    end

    # Optional Selective Dynamics
    if lowercase(lines[idx]) |> x -> startswith(x,"s")
        idx += 1
    end

    coord_type = lowercase(lines[idx])
    idx += 1

    natoms = sum(counts)

    coords = zeros(natoms,3)

    for i in 1:natoms
        fields = split(lines[idx])

        coords[i,:] .= parse.(Float64, fields[1:3])

        idx += 1
    end

    return POSCAR(
        lattice,
        symbols,
        coords,
        coord_type
    )
end

