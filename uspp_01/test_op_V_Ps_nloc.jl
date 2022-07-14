using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../pwscf_02/PWSCFInput.jl")

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
    
    ik = 1
    Ham.ik = ik
    Ngwk = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates

    psi = zeros(ComplexF64, Ngwk, Nstates)
    for i in 1:Nstates
        psi[i,i] = 1.0
    end

    println("sum psi = ", sum(psi))

    Vpsi = zeros(ComplexF64, Ngwk, Nstates)
    op_V_Ps_nloc!(Ham, psi, Vpsi)

    println("Sum Vpsi = ", sum(Vpsi))

end

test_main()
