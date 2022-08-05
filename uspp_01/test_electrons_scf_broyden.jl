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

include("electrons_scf_broyden.jl")
include("mix_broyden_02.jl")
#include("mix_broyden_03.jl")

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



function prepare_scf!(Ham)
    # Initial density
    Rhoe, RhoeG = atomic_rho_g(Ham)
    # Update the potentials
    _ = update_from_rhoe!( Ham, Rhoe, RhoeG )
    # Initial wavefunc
    psiks = rand_BlochWavefunc(Ham)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    # Reorthonormalize (FIXME: need to be included in rand_Blochwavefunc)
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin - 1)*Nkpt
        ortho_sqrt_with_S!(Ham, psiks[ikspin])
    end
    return psiks
end



function test_main()

    Ham = init_Ham_from_pwinput()
    println(Ham)

    psiks = prepare_scf!(Ham)

    electrons_scf_broyden!(Ham, psiks, NiterMax=100)

end

test_main()
