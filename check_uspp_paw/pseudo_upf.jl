function load_Dion( isp::Int64 )
    json_data = JSON.parsefile("pseudo_upf_"*string(isp)*".json")["pseudo_upf"];
    shape_dion = convert(Vector{Int64}, json_data["shape_dion"])
    dion = reshape(
        convert(Vector{Float64}, json_data["dion"]), shape_dion...
    )*2 # convert to 1/Ry to 1/Ha
    return dion
end

function load_qfuncl( isp::Int64 )
    json_data = JSON.parsefile("pseudo_upf_"*string(isp)*".json")["pseudo_upf"];
    shape_qfuncl = convert(Vector{Int64}, json_data["shape_qfuncl"])
    qfuncl = reshape(convert(Vector{Float64}, json_data["qfuncl"]), shape_qfuncl...)
    return qfuncl
end

function load_psp_vloc( isp::Int64 )
    json_data = JSON.parsefile("pseudo_upf_"*string(isp)*".json")["pseudo_upf"];
    vloc = convert(Vector{Float64}, json_data["vloc"])*0.5
    return vloc
end

function check_pseudo_upf(Ham)
    Nspecies = Ham.atoms.Nspecies

    for isp in 1:Nspecies

        println()
        println("---------------------------------")
        println("Checking pseudo_upf for ", Ham.pspots[isp].atsymb)
        println("---------------------------------")
        
        dion = load_Dion(isp)
        diff1 = Ham.pspots[isp].Dion - dion
        println("Dion: maximum abs diff = ", maximum(abs.(diff1)))

        vloc = load_psp_vloc(isp)
        diff1 = Ham.pspots[isp].V_local - vloc
        println("V_local: maximum abs diff = ", maximum(abs.(diff1)))

        if Ham.pspots[isp].is_ultrasoft
            qfuncl = load_qfuncl(isp)
            diff1 = Ham.pspots[isp].qfuncl - qfuncl
            println("qfuncl: maximum abs diff = ", maximum(abs.(diff1)))
        end

    end

end


