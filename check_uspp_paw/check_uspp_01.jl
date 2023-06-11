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
    # Initial density
    Rhoe, RhoeG = atomic_rho_g(Ham);
    # Also initialize becsum in case of PAW
    if any(Ham.pspotNL.are_paw)
        PAW_atomic_becsum!(Ham);
    end
    # Update the potentials
    Ehartree, Exc, Evtxc = update_from_rhoe!( Ham, nothing, Rhoe, RhoeG )
    Ham.energies.Hartree = Ehartree
    Ham.energies.XC = Exc
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

function load_qq_at_nt()
    json_data = JSON.parsefile("uspp_mod.json")["uspp"]
    
    shape_qq_nt = convert(Vector{Int64}, json_data["shape_qq_nt"])
    qq_nt = reshape(
        convert(Vector{Float64}, json_data["qq_nt"]),
        shape_qq_nt...
    )

    shape_qq_at = convert(Vector{Int64}, json_data["shape_qq_at"])
    qq_at = reshape(
        convert(Vector{Float64}, json_data["qq_at"]),
        shape_qq_at...
    )

    return qq_at, qq_nt
end



function check_qq_at_qq_nt(Ham::Hamiltonian)

    if !any(Ham.pspotNL.are_ultrasoft)
        println("No USPP or PAW data: skipped")
        return
    end

    println()
    println("------------------------")
    println("Checking qq_at and qq_nt")
    println("------------------------")

    Q_CONV_FACTOR = 4

    qq_at, qq_nt = load_qq_at_nt()

    println()
    println("qq_at")
    diff1 = Ham.pspotNL.qq_at - qq_at*Q_CONV_FACTOR
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))

    println()
    println("qq_nt")
    diff1 = Ham.pspotNL.qq_nt - qq_nt*Q_CONV_FACTOR
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))
end


function load_dvan()
    json_data = JSON.parsefile("uspp_mod.json")["uspp"]
    shape_dvan = convert(Vector{Int64}, json_data["shape_dvan"])
    dvan = reshape(
        convert(Vector{Float64}, json_data["dvan"]),
        shape_dvan...
    )
    return dvan
end

function check_Dvan(Ham::Hamiltonian)

    if !any(Ham.pspotNL.are_ultrasoft)
        println("No USPP or PAW data: skipped")
        return
    end

    println()
    println("-------------")
    println("Checking Dvan")
    println("-------------")

    CONV_FACTOR = 2

    dvan = load_dvan()

    diff1 = Ham.pspotNL.Dvan - dvan*CONV_FACTOR
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))
end


function subtract_Dvan(Ham::Hamiltonian)

    Deeq = Ham.pspotNL.Deeq
    Dvan = Ham.pspotNL.Dvan

    Deeq_new = copy(Deeq)
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    Nspin = size(Deeq, 4)
    nh = Ham.pspotNL.nh

    for ia in 1:Natoms
        isp = atm2species[ia]
        for ispin in 1:Nspin
            for ih in 1:nh[isp], jh in ih:nh[isp]
                Deeq_new[ih,jh,ia,ispin] = Deeq[ih,jh,ia,ispin] - Dvan[ih,jh,isp]
                Deeq_new[jh,ih,ia,ispin] = Deeq_new[ih,jh,ia,ispin]
            end
        end
    end
    return Deeq_new
end

# Deeq and Dvan can be from QE
function subtract_Dvan(Ham::Hamiltonian, Deeq, Dvan)

    Deeq_new = copy(Deeq)
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    Nspin = size(Deeq, 4)
    nh = Ham.pspotNL.nh

    for ia in 1:Natoms
        isp = atm2species[ia]
        for ispin in 1:Nspin
            for ih in 1:nh[isp], jh in ih:nh[isp]
                Deeq_new[ih,jh,ia,ispin] = Deeq[ih,jh,ia,ispin] - Dvan[ih,jh,isp]
                Deeq_new[jh,ih,ia,ispin] = Deeq_new[ih,jh,ia,ispin]
            end
        end
    end
    return Deeq_new
end



function load_deeq()
    json_data = JSON.parsefile("uspp_mod.json")["uspp"]
    shape_deeq = convert(Vector{Int64}, json_data["shape_deeq"])
    deeq = reshape(
        convert(Vector{Float64}, json_data["deeq"]),
        shape_deeq...
    )
    return deeq
end


function check_Deeq(Ham::Hamiltonian)

    if !any(Ham.pspotNL.are_ultrasoft)
        println("No USPP or PAW data: skipped")
        return
    end

    println()
    println("-------------")
    println("Checking Deeq")
    println("-------------")

    CONV_FACTOR = 2

    deeq = load_deeq()

    diff1 = Ham.pspotNL.Deeq - deeq*CONV_FACTOR
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))
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
    vltot = convert(Vector{Float64}, json_data["vltot"])
    return vltot*0.5 # convert to Ha
end

function check_V_Ps_loc(Ham)
    println()
    println("-----------------")
    println("Checking V_Ps_loc")
    println("-----------------")

    vltot = load_V_Ps_loc()
    diff1 = Ham.potentials.Ps_loc - vltot
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


function load_vrs()
    json_data = JSON.parsefile("scf_mod.json")["scf_mod"]
    shape_vrs = convert(Vector{Int64}, json_data["shape_vrs"])
    vrs = reshape(
        convert(Vector{Float64}, json_data["vrs"]),
        shape_vrs...
    )
    return vrs*0.5 # convert to Ha
end

function load_v_of_r()
    json_data = JSON.parsefile("scf_mod.json")["scf_mod"]
    shape_v_of_r = convert(Vector{Int64}, json_data["shape_v_of_r"])
    v_of_r = reshape(
        convert(Vector{Float64}, json_data["v_of_r"]),
        shape_v_of_r...
    )
    return v_of_r*0.5 # convert to Ha
end


function check_v_of_r(Ham)

    println()
    println("---------------")
    println("Checking v_of_r")
    println("---------------")

    v_of_r = load_v_of_r()
    diff1 = (Ham.potentials.XC + Ham.potentials.Hartree) - v_of_r
    println("maximum abs diff = ", maximum(abs.(diff1)))
    println("minimum abs diff = ", minimum(abs.(diff1)))
    println("mean    abs diff = ", mean(abs.(diff1)))
end


function load_rho_of_r()
    json_data = JSON.parsefile("scf_mod.json")["scf_mod"]
    shape_rho_of_r = convert(Vector{Int64}, json_data["shape_rho_of_r"])
    rho_of_r = reshape(
        convert(Vector{Float64}, json_data["rho_of_r"]),
        shape_rho_of_r...
    )
    return rho_of_r
end


function check_Rhoe(Ham)
    println()
    println("-------------")
    println("Checking Rhoe")
    println("-------------")

    rho_of_r = load_rho_of_r()
    diff1 = Ham.rhoe - rho_of_r
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
    check_Rhoe(Ham)
    check_v_of_r(Ham)
    check_clebsch_gordan(Ham)

    return
end

main()