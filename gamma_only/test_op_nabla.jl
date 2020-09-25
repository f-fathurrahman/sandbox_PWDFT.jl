using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("ortho_GS_gamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")
include("BlochWavefuncGamma.jl")
include("calc_rhoe_gamma.jl")
include("Poisson_solve_gamma.jl")
include("op_K_gamma.jl")
include("op_V_loc_gamma.jl")
include("op_V_Ps_nloc_gamma.jl")
include("op_H_gamma.jl")

include("unfold_BlochWavefuncGamma.jl")

function test_op_nabla()

    Random.seed!(1234)

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
    #               LatVecs = gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "NH3.xyz") )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    # Initialize Hamiltonian
    ecutwfc = 50.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )

    # Compare the size of Hamiltonian
    println("Size Ham gamma   = ", Base.summarysize(Ham) /1024/1024)
    println("Size Ham (usual) = ", Base.summarysize(Ham_)/1024/1024)

    psis = randn_BlochWavefuncGamma(Ham)
    ortho_check(psis)

    psis2 = randn_BlochWavefuncGamma(Ham)
    ortho_check(psis2)

    Nstates = size(psis.data[1],2)
    ispin = 1
    println("Some psis.data")
    for ist in 1:Nstates
        @printf("State #%d\n", ist)
        for igw in 1:4
            c = psis.data[ispin][igw,ist]
            @printf("%3d %3d [%18.10f,%18.10f]\n", igw, ist, c.re, c.im)
        end
    end

    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
    psiks2 = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis2 )

    println("\nSome psiks[1]")
    for ist in 1:Nstates
        @printf("State #%d\n", ist)
        for igw in 1:4
            c = psiks[ispin][igw,ist]
            @printf("%3d %3d [%18.10f,%18.10f]\n", igw, ist, c.re, c.im)
        end
    end

    println("\ndot the same wavefunc")
    #
    @time s1 = dot(psiks,psiks)
    @time s1 = dot(psiks,psiks)
    println("Using dot (usual)           : ", s1)
    #
    @time s2 = dot_BlochWavefuncGamma(psis,psis)
    @time s2 = dot_BlochWavefuncGamma(psis,psis)
    println("Using dot_BlochWavefuncGamma: ", s2)
    println("diff: ", s1-s2, " (should be small)")

    println("\ndot the different wavefunc")
    #
    @time s1 = dot(psiks,psiks2)
    @time s1 = dot(psiks,psiks2)
    println("Using dot (usual)           : ", s1)
    #
    @time s2 = dot_BlochWavefuncGamma(psis,psis2)
    @time s2 = dot_BlochWavefuncGamma(psis,psis2)
    println("Using dot_BlochWavefuncGamma: ", s2)
    println("diff: ", s1-s2, " (should be small)")

    println("\nCompare calc_rhoe:")
    #
    println("Using Gamma-point trick")
    Rhoe = calc_rhoe(Ham, psis); print("time calc_rhoe: ")
    @time Rhoe = calc_rhoe(Ham, psis)
    #
    println("Using usual kpt")
    Rhoe_ = calc_rhoe(Ham_, psiks); print("time calc_rhoe: ")
    @time Rhoe = calc_rhoe(Ham_, psiks)


    Npoints = prod(Ham.pw.Ns)
    println("\nCompare op_nabla:")
    #
    println("Using Gamma-point trick")
    RhoeTot = dropdims(sum(Rhoe, dims=2),dims=2)
    ∇Rhoe = op_nabla(Ham.pw, RhoeTot); print("time op_nabla: ")
    @time ∇Rhoe = op_nabla(Ham.pw, RhoeTot)
    #
    println("Using usual kpt")
    RhoeTot_ = dropdims(sum(Rhoe_, dims=2),dims=2)
    ∇Rhoe_ = op_nabla(Ham_.pw, RhoeTot_); print("time op_nabla: ")
    @time ∇Rhoe_ = op_nabla(Ham_.pw, RhoeTot_)

    s1 = sum(∇Rhoe)
    s2 = sum(∇Rhoe_)
    println("s1 = ", s1)
    println("s2 = ", s2)
    #
    mae = sum(abs.(∇Rhoe .- ∇Rhoe_))/Npoints
    println("diff: ", mae, " (should be small)")

end

test_op_nabla()