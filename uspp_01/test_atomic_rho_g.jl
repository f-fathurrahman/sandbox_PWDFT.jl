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
    if abs(dual-4) >= 1e-10
        error("Please use NCPP only")
    end

    pspfiles = pwinput.pspfiles
    # Need special treatement for GTH ?

    return Hamiltonian(atoms, pspfiles, ecutwfc)
end


function test_main()
    
    Ham = init_Ham_from_pwinput()
    
    Rhoe = atomic_rho_g(Ham)

end

test_main()
