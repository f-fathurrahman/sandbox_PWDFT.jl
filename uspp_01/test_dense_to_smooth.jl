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

include("dense_to_smooth.jl")

function test_main()
    
    Ham = init_Ham_from_pwinput()

    @assert Ham.pw.using_dual_grid

    NptsDense = prod(Ham.pw.Ns)
    NptsSmooth = prod(Ham.pw.Nss)

    vin = zeros(Float64, NptsDense)
    vout = zeros(Float64, NptsSmooth)

    vin[:] .= 1.1
    vin[1:10] .= 2.5

    println("sum vin before dense_to_smooth = ", sum(vin))
    println("sum vout before dense_to_smooth = ", sum(vout))

    dense_to_smooth!( Ham.pw, vin, vout )

    println("sum vin after dense_to_smooth = ", sum(vin))
    println("sum vout after dense_to_smooth = ", sum(vout))

    println("Some vin and vout")
    for i in 1:10
        @printf("%5d %18.10f %18.10f\n", i, vin[i], vout[i])
    end

end

test_main()
