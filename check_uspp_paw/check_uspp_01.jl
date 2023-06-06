# Need data from export2pwdftjl

using Printf
using OffsetArrays
using LinearAlgebra
using Serialization
import Statistics: mean

using Random
Random.seed!(1234)

using JSON

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

function create_Ham(filename)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename);
    return Ham
end

function load_clebsch_gordan()
    json_data = JSON.parsefile("uspp_mod.json")["uspp"]
    
    shape_ap = convert(Vector{Int64}, json_data["shape_ap"])
    ap = reshape(
        convert(Vector{Float64}, json_data["ap"]),
        shape_ap...
    )

    shape_lpx = convert(Vector{Int64}, json_data["shape_lpx"])
    lpx = reshape(
        convert(Vector{Int64}, json_data["lpx"]),
        shape_lpx...
    )

    shape_lpl = convert(Vector{Int64}, json_data["shape_lpl"])
    lpl = reshape(
        convert(Vector{Int64}, json_data["lpl"]),
        shape_lpl...
    )

    return ap, lpx, lpl
end


function load_qradG()
    json_data = JSON.parsefile("uspp_mod.json")["uspp"]
    
    shape_qrad = convert(Vector{Int64}, json_data["shape_qrad"])
    qradG = reshape(
        convert(Vector{Float64}, json_data["qrad"]),
        shape_qrad...
    )
    return qradG
end




function check_clebsch_gordan(Ham)

    println()
    println("------------------------------")
    println("Checking Clebsch-Gordan arrays")
    println("------------------------------")

    if !any(Ham.pspotNL.are_ultrasoft)
        println("No USPP or PAW data: skipped")
        return
    end

    ap, lpx, lpl = load_clebsch_gordan()

    diff1 = abs.(Ham.pspotNL.ap - ap)
    println("diff abs ap: ", mean(diff1))

    diff1 = abs.(Ham.pspotNL.lpx - lpx)
    println("diff abs lpx: ", mean(diff1))

    diff1 = abs.(Ham.pspotNL.lpl - lpl)
    println("diff abs lpl: ", mean(diff1))

end




# FIXME: Not yet checked
function load_vkb()
    json_data = JSON.parsefile("uspp_mod.json")["uspp"]
    
    shape_vkb = convert(Vector{Int64}, json_data["shape_vkb"])
    
    vkb_real = reshape(
        convert(Vector{Float64}, json_data["vkb_real"]),
        shape_vkb...
    )

    vkb_imag = reshape(
        convert(Vector{Float64}, json_data["vkb_imag"]),
        shape_vkb...
    )

    vkb = vkb_real + im*vkb_imag

    return vkb
end


function load_V_Ps_loc()
    json_data = JSON.parsefile("scf_mod.json")["scf_mod"]
    vltot = convert(Vector{Float64}, json_data["vltot"])*0.5 # convert to Ha
    return vltot
end

function check_V_Ps_loc(Ham)
    println()
    println("-----------------")
    println("Checking V_Ps_loc")
    println("-----------------")

    vltot_Ha = load_V_Ps_loc()
    diff1 = Ham.potentials.Ps_loc  - vltot_Ha
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))
end


function load_rhoe_core()
    json_data = JSON.parsefile("scf_mod.json")["scf_mod"]
    rhoe_core = convert(Vector{Float64}, json_data["rho_core"])
    return rhoe_core
end


function check_rhoe_core(Ham)

    println()
    println("------------------")
    println("Checking rhoe_core")
    println("------------------")

    Nspecies = Ham.atoms.Nspecies
    pspots = Ham.pspots
    if eltype(Ham.pspots) == PsPot_GTH
        println("Checking rhoe_core is skipped: using PsPot_GTH")
        return
    end

    are_psp_using_nlcc = zeros(Bool,Nspecies)
    for isp in 1:Nspecies
        if pspots[isp].is_nlcc
            are_psp_using_nlcc[isp] = true
        end
    end
    if !any(are_psp_using_nlcc)
        println("Checking rhoe_core is skipped: no core-correction")
        return
    end

    rhoe_core = load_rhoe_core()
    diff1 = Ham.rhoe_core  - rhoe_core
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))

    CellVolume = Ham.pw.CellVolume
    Npoints = prod(Ham.pw.Ns)
    rhoneg = 0.0
    for x in rhoe_core
        if x < 0.0
            rhoneg += abs(x)
        end
    end
    println("Check negative rhoe_core = ", rhoneg*CellVolume/Npoints)
    println("integ rhoe_core = ", sum(rhoe_core)*CellVolume/Npoints)

end





function load_Vloc_gl_sp()
    json_data = JSON.parsefile("vlocal_mod.json")["vlocal"]
    shape_vloc = convert(Vector{Int64}, json_data["shape_vloc"])
    Vloc_gl_sp = reshape(
        convert(Vector{Float64}, json_data["vloc"]),
        shape_vloc...
    )
    return Vloc_gl_sp*0.5 # convert to Ha
end


function check_Vloc_gl_sp(Ham)

    println()
    println("-------------------")
    println("Checking Vloc_gl_sp")
    println("-------------------")

    # This quantity is not stored. We recalculate it here
    pspots = Ham.pspots
    Nspecies = Ham.atoms.Nspecies
    pw = Ham.pw
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl_sp = zeros(Float64, Ngl, Nspecies)
    for isp in 1:Nspecies
        @views eval_Vloc_G!( pspots[isp], G2_shells, Vgl_sp[:,isp] )
    end

    # Scale by pw.CellVolume
    Vgl_sp[:,:] /= pw.CellVolume

    Vloc_gl_sp = load_Vloc_gl_sp()

    diff1 = Vgl_sp - Vloc_gl_sp
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))

end


include("pseudo_upf.jl")

function main()
    if length(ARGS) >= 1
        filename = ARGS[1]
    else
        filename = "PWINPUT"
    end
    
    Ham = create_Ham(filename)

    check_pseudo_upf(Ham)
    check_rhoe_core(Ham)
    check_V_Ps_loc(Ham)
    check_Vloc_gl_sp(Ham)
    check_clebsch_gordan(Ham)


    return
end

main()