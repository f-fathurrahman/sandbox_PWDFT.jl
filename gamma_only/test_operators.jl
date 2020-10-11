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

function compare_potentials(Ham, Ham_)

    println("\nV Ps loc comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.Ps_loc[ip], Ham_.potentials.Ps_loc[ip])
    end

    println("\nV Hartree comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.Hartree[ip], Ham_.potentials.Hartree[ip])
    end

    println("\nV XC comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.XC[ip,1], Ham_.potentials.XC[ip,1])
    end

    println("\nV Total comparison")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, Ham.potentials.Total[ip,1], Ham_.potentials.Total[ip,1])
    end

    return
end

function test_01()

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


    update!(Ham, Rhoe)
    update!(Ham_, Rhoe_)
    #compare_potentials(Ham, Ham_)

    #
    # Compare the result of applying operators
    #

    println("\nCompare op_K result")
    println("---------------------\n")
    #
    println("Using Gamma-point trick")
    Kpsis = op_K(Ham, psis); print("time op_K: ")
    @time Kpsis = op_K(Ham, psis)
    s1 = dot_BlochWavefuncGamma(psis, Kpsis)
    #
    println("Using usual kpt")
    Kpsiks = op_K(Ham_, psiks); print("time op_K: ")
    @time Kpsiks = op_K(Ham_, psiks)
    s2 = dot(psiks, Kpsiks)
    #
    println("dot(psis,Kpsis)   = ", s1)
    println("dot(psiks,Hpsiks) = ", s2)
    println("diff = ", s1 - s1, " (should be small)")

    println("Some Kpsi: (some symmetry should exist)")
    for ig in 1:5
        v1r = real(Kpsis.data[1][ig,1])
        v1i = imag(Kpsis.data[1][ig,1])
        v2r = real(Kpsiks[1][ig,1])
        v2i = imag(Kpsiks[1][ig,1])
        @printf("%5d (%18.10f,%18.10f) (%18.10f,%18.10f)\n", ig, v1r, v1i, v2r, v2i)
    end

    println("\nCompare op_V_loc")
    println("------------------\n")
    #
    println("Using Gamma-point trick")
    V_loc_psis = op_V_loc(Ham, psis); print("time op_V_loc: ")
    @time V_loc_psis = op_V_loc(Ham, psis)
    s1 = dot_BlochWavefuncGamma(psis, V_loc_psis)

    println("Using usual kpt")
    V_loc_psiks = op_V_loc(Ham_, psiks); print("time op_V_loc: ")
    @time V_loc_psiks = op_V_loc(Ham_, psiks)
    s2 = dot(psiks, V_loc_psiks)
    
    println("dot V_loc_psis  = ", s1)
    println("dot V_loc_psiks = ", s2)
    println("diff = ", s1 - s1, " (should be small)")
    
    println("Some Vloc_psi: (some symmetry should exist)")
    for ig in 1:5
        v1r = real(V_loc_psis.data[1][ig,1])
        v1i = imag(V_loc_psis.data[1][ig,1])
        v2r = real(V_loc_psiks[1][ig,1])
        v2i = imag(V_loc_psiks[1][ig,1])
        @printf("%5d (%18.10f,%18.10f) (%18.10f,%18.10f)\n", ig, v1r, v1i, v2r, v2i)
    end

    println("\nCompare op_V_Ps_nloc")
    println("----------------------\n")
    #
    println("Using Gamma-point trick")
    V_psis = op_V_Ps_nloc(Ham, psis); print("time op_V_Ps_nloc: ")
    @time V_psis = op_V_Ps_nloc(Ham, psis)
    s1 = dot_BlochWavefuncGamma(psis, V_psis)

    println("Using usual kpt")
    V_psiks = op_V_Ps_nloc(Ham_, psiks); print("time op_V_Ps_nloc: ")
    @time V_psiks = op_V_Ps_nloc(Ham_, psiks)
    s2 = dot(psiks, V_psiks)
    
    println("dot V_Ps_nloc_psis  = ", s1)
    println("dot V_Ps_nloc_psiks = ", s2)
    println("diff = ", s1 - s1, " (should be small)")

    println("Some V_Ps_nloc_psi: (some symmetry should exist)")
    for ig in 1:5
        v1r = real(V_psis.data[1][ig,1])
        v1i = imag(V_psis.data[1][ig,1])
        v2r = real(V_psiks[1][ig,1])
        v2i = imag(V_psiks[1][ig,1])
        @printf("%5d (%18.10f,%18.10f) (%18.10f,%18.10f)\n", ig, v1r, v1i, v2r, v2i)
    end

    println("\nCompare op_H")
    println("-------------\n")
    #
    println("Using Gamma-point trick")
    Hpsis = op_H(Ham, psis); print("time op_H: ")
    @time Hpsis = op_H(Ham, psis)
    s1 = dot_BlochWavefuncGamma(psis, Hpsis)

    println("Using usual kpt")
    Hpsiks = op_H(Ham_, psiks); print("time op_H: ")
    @time Hpsiks = op_H(Ham_, psiks)
    s2 = dot(psiks, Hpsiks)
    
    println("dot Hpsis  = ", s1)
    println("dot Hpsiks = ", s2)
    println("diff = ", s1 - s1, " (should be small)")

end

test_01()