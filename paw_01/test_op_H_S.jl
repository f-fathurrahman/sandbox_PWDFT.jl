using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

function test_main()

    Ham = init_Ham_from_pwinput()
    println(Ham)

    Rhoe, RhoeG = atomic_rho_g(Ham)
    if any(Ham.pspotNL.are_paw)
        PAW_atomic_becsum!(Ham)
    end

    Ehartree, Exc, Evtxc = update_from_rhoe!( Ham, nothing, Rhoe, RhoeG )

    ik = 1
    Ham.ik = ik # set the current k index of H
    Ngwk = Ham.pw.gvecw.Ngw[ik]
    Nstates = Ham.electrons.Nstates

    # Prepare wavefunctions (not necessarily normalized w.r.t op_S)
    psi = zeros(ComplexF64, Ngwk, Nstates)
    for i in 1:Nstates
        psi[i,i] = 1.1
    end

    Spsi = zeros(ComplexF64, Ngwk, Nstates)
    Hpsi = zeros(ComplexF64, Ngwk, Nstates)

    println("Before: sum(psi) = ", sum(psi))
    println("Before: sum(Hpsi) = ", sum(Hpsi))
    println("Before: sum(Spsi) = ", sum(Spsi))

    op_H!(Ham, psi, Hpsi)
    op_S!(Ham, psi, Spsi)

    println("\nTry 1")
    println("sum Hpsi = ", sum(Hpsi))
    println("sum(abs(Hpsi)) = ", sum(abs.(Hpsi)))
    println()
    println("sum Spsi = ", sum(Spsi))
    println("sum(abs(Spsi)) min ref = ", sum(abs.(Spsi)) - sum(psi))

    fill!(Hpsi, 0.0)
    op_H!(Ham, psi, Hpsi)

    println("\nTry 2")
    println("sum Hpsi = ", sum(Hpsi))
    println("sum(abs(Hpsi)) = ", sum(abs.(Hpsi)))

end

test_main()
