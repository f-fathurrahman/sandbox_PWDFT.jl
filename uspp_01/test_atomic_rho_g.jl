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


function test_main()
    
    Ham = init_Ham_from_pwinput()

    println("Npoints = ", prod(Ham.pw.Ns))
    println("Ns = ", Ham.pw.Ns)
    
    if Ham.pw.using_dual_grid
        println("Npoints smooth = ", prod(Ham.pw.Nss))
        println("Nss = ", Ham.pw.Nss)
    end

    Rhoe, RhoeG = atomic_rho_g(Ham)

end

test_main()
