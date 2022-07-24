using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("update_from_rhoe.jl")
include("newd.jl")
include("op_S.jl")

#atoms, pw, pspots = init_test_main()
function test_main()

    Ham = init_Ham_from_pwinput()


    Rhoe, RhoeG = atomic_rho_g(Ham)

    update_from_rhoe!( Ham, Rhoe, RhoeG )

    calc_newDeeq!( Ham )

    ik = 1
    Ham.ik = ik # set the current k index of H
    Ngwk = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates

    psi = zeros(ComplexF64, Ngwk, Nstates)
    for i in 1:Nstates
        psi[i,i] = 1.0
    end

    Spsi = zeros(ComplexF64, Ngwk, Nstates)
    Hpsi = zeros(ComplexF64, Ngwk, Nstates)

    op_H!(Ham, psi, Hpsi)
    op_S!(Ham, psi, Spsi)

    println("sum Hpsi = ", sum(Hpsi))
    println("sum(abs(Hpsi)) = ", sum(abs.(Hpsi)))

    println("sum Spsi = ", sum(Spsi))
    println("sum(abs(Spsi)) min ref = ", sum(abs.(Spsi)) - sum(psi))

end

test_main()
