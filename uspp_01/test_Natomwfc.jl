using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("newd.jl")

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
    for isp in 1:atoms.Nspecies
        if is_using_extension_gth(pspfiles[isp])
            error("GTH pspot is not yet supported")
        end
    end

    return Hamiltonian(atoms, pspfiles, ecutwfc, dual=dual)
end



function test_main()
    
    Ham = init_Ham_from_pwinput()
    Natomwfc = calc_Natomwfc(Ham.atoms, Ham.pspots)
    println("Natomwfc = ", Natomwfc)

end

test_main()