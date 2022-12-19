using Printf
using OffsetArrays
using LinearAlgebra
using SpecialFunctions: sphericalbesselj

using Random
Random.seed!(1234)

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("smooth_to_dense.jl")
include("update_from_rhoe.jl")
include("newd.jl")
include("op_S.jl")
include("calc_rhoe_uspp.jl")
include("../diag_davidson_qe/diag_davidson_qe_v2.jl")

#include("electrons_scf_broyden.jl")
include("electrons_scf_broyden_v2.jl")
include("mix_broyden_02.jl")
#include("mix_broyden_03.jl")

include("ortho_with_S.jl")

function test_main()
    Ham = init_Ham_from_pwinput()
    write_xsf("ATOMS_from_pwinput.xsf", Ham.atoms)
    println(Ham)
    psiks = rand_BlochWavefunc(Ham)
    electrons_scf_broyden!(Ham, psiks, NiterMax=100)
end

test_main()
