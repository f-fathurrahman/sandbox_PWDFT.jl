function load_ddd_paw()
    json_data = JSON.parsefile("paw_variables.json")["paw_variables"]
    
    shape_ddd_paw = convert(Vector{Int64}, json_data["shape_ddd_paw"])
    ddd_paw = reshape(
        convert(Vector{Float64}, json_data["ddd_paw"]),
        shape_ddd_paw...
    )
    return ddd_paw*0.5 # convert to Ha
end


function check_ddd_paw(Ham)

    println()
    println("----------------")
    println("Checking ddd_paw")
    println("----------------")

    if !any(Ham.pspotNL.are_paw)
        println("No PAW data: skipped")
        return
    end

    ddd_paw = load_ddd_paw()

    diff1 = Ham.pspotNL.paw.ddd_paw - ddd_paw
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))
end


function load_paw_radial_scalars( isp::Int64 )
    json_data = JSON.parsefile("paw_radial_integrator_"*string(isp)*".json")["paw_radial_integrator"]
    lmax = json_data["lmax"]
    lm_max = json_data["lm_max"]
    ladd = json_data["ladd"]
    nx = json_data["nx"]
    return lmax, lm_max, ladd, nx
end

function load_paw_radial_ww( isp::Int64 )
    json_data = JSON.parsefile("paw_radial_integrator_"*string(isp)*".json")["paw_radial_integrator"]
    ww = json_data["ww"]
    return ww
end

function load_paw_radial_wwylm( isp::Int64 )
    json_data = JSON.parsefile("paw_radial_integrator_"*string(isp)*".json")["paw_radial_integrator"]
    shape_wwylm = convert(Vector{Int64}, json_data["shape_wwylm"])
    wwylm = reshape(
        convert(Vector{Float64}, json_data["wwylm"]),
        shape_wwylm...
    )
    return wwylm
end

function check_paw_radial(Ham)

    Nspecies = Ham.atoms.Nspecies

    for isp in 1:Nspecies

        println()
        println("-------------------------------------------------")
        println("Checking PAW radial integator for species ", isp)
        println("-------------------------------------------------")

        if !Ham.pspotNL.are_paw[isp]
            println("No PAW data: skipped")
            continue
        end

        # Access PAWAtomicSphere for this species
        paw_r = Ham.pspotNL.paw.spheres[isp]

        lmax, lm_max, ladd, nx = load_paw_radial_scalars(isp)
        
        diff1 = paw_r.lmax - lmax
        println("lmax: diff = ", diff1)
        
        diff1 = paw_r.lm_max - lm_max
        println("lm_max: diff = ", diff1)
        
        diff1 = paw_r.ladd - ladd
        println("ladd: diff = ", diff1)
        
        diff1 = paw_r.nx - nx
        println("lmax: diff = ", diff1)

        ww = load_paw_radial_ww(isp)
        diff1 = paw_r.ww - ww
        println("ww: maximum abs diff = ", maximum(abs.(diff1)))

        wwylm = load_paw_radial_wwylm(isp)
        diff1 = paw_r.wwylm - wwylm
        println("ww: maximum abs diff = ", maximum(abs.(diff1)))

    end

end


