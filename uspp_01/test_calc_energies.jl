using Printf
using OffsetArrays
using LinearAlgebra
using SpecialFunctions: sphericalbesselj
using Random

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("update_from_rhoe.jl")
include("newd.jl")
include("op_S.jl")
include("my_scf_01.jl")
include("calc_rhoe_uspp.jl")
include("../diag_davidson_qe/diag_davidson_qe_v2.jl")

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
    for ist = 2:Nstates
        c = dot( psi[:,ist], Spsi[:,1] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end


Random.seed!(1234)

function test_main()

    Ham = init_Ham_from_pwinput()

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    # Initial density
    Rhoe, RhoeG = atomic_rho_g(Ham)

    Ehartree, Exc, Evtxc = update_from_rhoe!( Ham, Rhoe, RhoeG )

    # Initial wavefunc
    psiks = rand_BlochWavefunc(Ham)
    # Reorthonormalize (FIXME: need to be included in rand_Blochwavefunc)
    for i in 1:size(psiks,1)
        ortho_sqrt_with_S!(Ham, psiks[i])
    end

    psi = psiks[1]

    Hpsi = op_H(Ham, psi)
    Spsi = op_S(Ham, psi)
    Kpsi = op_K(Ham, psi)
    Vnlpsi = op_V_Ps_nloc(Ham, psi)

    energies = calc_energies(Ham, psiks)

    println(energies)

    println("dot psi Hpsi = ", dot(psi, Hpsi))
    println("dot psi Spsi = ", dot(psi, Spsi))

    println("dot psi Kpsi = ", 2*dot(psi, Kpsi)/dot(psi,psi))

    @assert Ham.pw.gvecw.kpoints.Nkpt == 1
    @assert Ham.electrons.Nspin == 1

    ikspin = 1
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Ekin = 0.0
    EpsNL = 0.0
    for ist in 1:Nstates
        Ekin += Focc[ist,ikspin] * dot(psi[:,ist], Kpsi[:,ist])
    end
    println("Ekin = ", Ekin)


    println("dot psi Vnlpsi = ", 2*dot(psi, Vnlpsi)/dot(psi,psi))

end

test_main()
