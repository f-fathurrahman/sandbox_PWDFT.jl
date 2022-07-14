using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("atomic_rho_g.jl")

function init_Ham_from_pwinput()
    println("ARGS = ", ARGS)
    @assert length(ARGS) == 1
    pwinput = PWSCFInput(ARGS[1])

    atoms = pwinput.atoms

    ecutwfc = pwinput.ecutwfc
    ecutrho = pwinput.ecutrho
    dual = ecutrho/ecutwfc

    pspfiles = pwinput.pspfiles
    # Need special treatement for GTH ?

    return Hamiltonian(atoms, pspfiles, ecutwfc, dual=dual)
end


function add_V_xc!(Ham, Rhoe, RhoeG, Veff)

    Nspin = Ham.electrons.Nspin
    
    @assert Nspin == 1

    Npoints = size(Rhoe, 1)
    Vxc = zeros(Float64, Npoints)
    epsxc = zeros(Float64, Npoints)

    # XC potential
    # VWN is the default
    if Ham.rhoe_core == nothing
        epsxc[:], Vxc[:] = calc_epsxc_Vxc_VWN( Ham.xc_calc, Rhoe[:,1] )
    else
        println("sum rhoe_core = ", sum(Ham.rhoe_core))
        println("sum rhoe + rhoe_core = ", sum(Rhoe[:,1] + Ham.rhoe_core))
        epsxc[:], Vxc[:] = calc_epsxc_Vxc_VWN( Ham.xc_calc, Rhoe[:,1] + Ham.rhoe_core )
    end

    println("sum abs epsxc = ", sum(abs.(epsxc)))
    println("sum abs Vxc   = ", sum(abs.(Vxc)))

    dVol = Ham.pw.CellVolume / Npoints
    # Calculate etxc and vtxc
    Exc = sum(epsxc .* (Rhoe[:,1] + Ham.rhoe_core))*dVol
    Vtxc = sum(Vxc .* Rhoe)*dVol # Vtxc does not include rhoe_core

    println("Exc  = ", Exc)
    println("Vtxc = ", Vtxc)

    return Exc, Vtxc
end


function test_main()
    
    Ham = init_Ham_from_pwinput()

    println("Npoints = ", prod(Ham.pw.Ns))
    println("Ns = ", Ham.pw.Ns)
    
    if Ham.pw.using_dual_grid
        println("Npoints smooth = ", prod(Ham.pw.Nss))
        println("Nss = ", Ham.pw.Nss)
    end

    Rhoe, RhoeG = atomic_rho_g(Ham)

    Nspin = Ham.electrons.Nspin
    Npoints = prod(Ham.pw.Ns)
    Veff = zeros(Float64, Npoints, Nspin)

    Exc, Vtxc = add_V_xc!( Ham, Rhoe, RhoeG, Veff )

end

test_main()
