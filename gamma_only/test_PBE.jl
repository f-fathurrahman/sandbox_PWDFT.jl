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


import PWDFT: calc_Vxc_PBE!
function calc_Vxc_PBE!(
    xc_calc::XCCalculator,
    pw::PWGridGamma,
    Rhoe::Array{Float64,1},
    Vxc
)

    Npoints = size(Rhoe,1)

    # calculate gRhoe2
    gRhoe = op_nabla(pw, Rhoe)
    gRhoe2 = zeros(Float64, Npoints)
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)

    for ip in 1:Npoints
        #
        ρ = Rhoe[ip]
        #
        _, vx = PWDFT.XC_x_slater( ρ )
        _, v1x, v2x = PWDFT.XC_x_pbe( ρ, gRhoe2[ip] )
        
        _, vc = PWDFT.XC_c_pw( ρ )
        _, v1c, v2c = PWDFT.XC_c_pbe( ρ, gRhoe2[ip] )

        Vxc[ip] = vx + vc + v1x + v1c
        
        v2xc = v2x + v2c
        hx[ip] = v2xc*gRhoe[1,ip]
        hy[ip] = v2xc*gRhoe[2,ip]
        hz[ip] = v2xc*gRhoe[3,ip]
    end

    dh = op_nabla_dot(pw, hx, hy, hz)

    for ip in 1:Npoints
        Vxc[ip] = Vxc[ip] - dh[ip]
    end

    return
end



function test_PBE()

    Random.seed!(1234)
    
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
    @time s2 = dot(psis,psis)
    @time s2 = dot(psis,psis)
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


    xc_calc = XCCalculator()
    Npoints = prod(Ham.pw.Ns)

    RhoeTot = dropdims(sum(Rhoe, dims=2),dims=2)
    RhoeTot_ = dropdims(sum(Rhoe_, dims=2),dims=2)

    Vxc = zeros(Float64,Npoints)
    Vxc_ = zeros(Float64,Npoints)

    println("\nCompare calc_Vxc_PBE!:")
    #
    println("Using Gamma-point trick")
    calc_Vxc_PBE!(xc_calc, Ham.pw, RhoeTot, Vxc); print("time calc_Vxc_PBE!: ")
    @time calc_Vxc_PBE!(xc_calc, Ham.pw, RhoeTot, Vxc)
    #
    println("Using usual kpt")
    calc_Vxc_PBE!(xc_calc, Ham_.pw, RhoeTot_, Vxc_); print("time calc_Vxc_PBE!: ")
    @time calc_Vxc_PBE!(xc_calc, Ham_.pw, RhoeTot_, Vxc_)

    println("sum Vxc  = ", sum(Vxc))
    println("sum Vxc_ = ", sum(Vxc_))
    mae = sum(abs.(Vxc .- Vxc_))/Npoints
    println("diff = ", mae, " (should be small)")

    println(Vxc[1:5])
    println(Vxc_[1:5])

end

test_PBE()