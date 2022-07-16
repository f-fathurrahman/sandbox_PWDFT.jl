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

    pspfiles = pwinput.pspfiles
    # Need special treatement for GTH ?

    return Hamiltonian(atoms, pspfiles, ecutwfc, dual=dual)
end

include("fft_interpolate.jl")

function test_main()
    
    Ham = init_Ham_from_pwinput()

    @assert Ham.pw.using_dual_grid

    NptsDense = prod(Ham.pw.Ns)

    vin = zeros(ComplexF64, NptsDense)

    vin[:] .= 1.0
    vin[1:10] .= 2.5

    println("sum vin = ", sum(vin))

    vout = copy(vin)
    #R_to_G!(Ham.pw, vout)
    #vout[:] /= NptsDense # scale to match QE convention
    #
    G_to_R!(Ham.pw, vout)
    vout[:] *= NptsDense # scale to match QE convention

    println("sum vout = ", sum(vout))

    println("Some vout")
    for i in 1:10
        @printf("%5d %18.10f %18.10f\n", i, vout[i].re, vout[i].im)
    end

end

test_main()
