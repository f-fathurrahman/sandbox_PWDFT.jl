struct ElkInput
    LatVecs::Matrix{Float64}
    species_files::Vector{String}
    natoms_per_species::Vector{Int64}
    atomic_positions::Vector{Matrix{Float64}}
end

function read_elk_input()
    f = open("elk.in", "r")
    lines = readlines(f)
    close(f)

    LatVecs = zeros(Float64, 3, 3)
    species_files = Vector{String}()
    natoms_per_species = Vector{Int64}()
    atomic_positions = Vector{Matrix{Float64}}()
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
                for ia in natmsp
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
    end

    return ElkInput(
        LatVecs, species_files,
        natoms_per_species, atomic_positions
    )
end
