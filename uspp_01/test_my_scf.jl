using Printf
using OffsetArrays
using LinearAlgebra
using SpecialFunctions: sphericalbesselj

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

include("my_scf_01.jl")
#include("my_scf_02.jl")

function ortho_sqrt_with_S!( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    O = psi' * op_S(Ham, psi)
    Udagger = inv(sqrt(O))
    psi[:,:] = psi*Udagger
    return
end


function ortho_check_with_S( Ham::Hamiltonian, psi::Array{ComplexF64,2}; dVol=1.0 )
    Nstates = size(psi,2)
    Spsi = op_S(Ham, psi)
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], Spsi[:,ist] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], Spsi[:,2] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end


function test_main()

    Ham = init_Ham_from_pwinput()

    println(Ham)

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    # Initial density
    Rhoe, RhoeG = atomic_rho_g(Ham)

    #println("integ Rhoe after atomic_rho_g (1): ", sum(Rhoe)*dVol)
    Ehartree, Exc, Evtxc = update_from_rhoe!( Ham, Rhoe, RhoeG )

    #println("integ Rhoe after atomic_rho_g (2): ", sum(Rhoe)*dVol)
    #println("After update_from_rhoe: ", sum(Ham.rhoe)*dVol)

    # Initial wavefunc
    psiks = rand_BlochWavefunc(Ham)
    # Reorthonormalize (FIXME: need to be included in rand_Blochwavefunc)
    for i in 1:size(psiks,1)
        ortho_sqrt_with_S!(Ham, psiks[i])
    end

    # Check
    ortho_check_with_S(Ham, psiks[1])
    #O = psiks[1]' * op_S(Ham, psiks[1])
    #println("Real O = ")
    #display(real.(O)); println()
    #println("Imag O = ")
    #display(imag.(O)); println()

    my_scf!(Ham, psiks, NiterMax=100)

end

test_main()
