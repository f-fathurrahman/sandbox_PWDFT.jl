using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_energies_grad.jl")
include("calc_energies_grad_gamma.jl")

function test_01()

    Random.seed!(1234)

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
    #               LatVecs = gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "NH3.xyz") )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )

    psis = randn_BlochWavefuncGamma(Ham)
    ortho_check(psis)

    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    Rhoe = calc_rhoe(Ham, psis)

    Rhoe_ = calc_rhoe(Ham_, psiks)

    update!(Ham, Rhoe)
    update!(Ham_, Rhoe_)

    ispin = 1

    Nstates = Ham.electrons.Nstates
    g = zeros(ComplexF64, Ham.pw.gvecw.Ngw, Nstates)
    g_ = zeros(ComplexF64, Ham_.pw.gvecw.Ngw[1], Nstates)

    Hsub = zeros(ComplexF64, Nstates, Nstates)
    Hsub_ = zeros(ComplexF64, Nstates, Nstates)

    calc_grad!(Ham, psis.data[1], g, Hsub)    
    print("Gamma-trick: ")
    @time calc_grad!(Ham, psis.data[1], g, Hsub)

    calc_grad!(Ham_, psiks[1], g_, Hsub_)
    print("Usual: ")    
    @time calc_grad!(Ham_, psiks[1], g_, Hsub_)

    display(Hsub); println()
    display(Hsub_); println()

    println("\nFirst component:")
    println("g[1,1]  = ", g[1,1])
    println("g_[1,1] = ", g_[1,1])

    println("\nTest dot")
    println("Using dot gamma  = ", dot_gamma(g,g))
    println("Using dot(g_,g_) = ", dot(g_,g_))
end

test_01()