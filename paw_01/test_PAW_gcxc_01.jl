using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

function main(;filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    becsum = PAW_atomic_becsum(Ham.atoms, Ham.pspots, Ham.pspotNL, Nspin=1)
    
    println("sum becsum before PAW_symmetrize: ", sum(becsum))
    PAW_symmetrize!(Ham, becsum)
    println("sum becsum after PAW_symmetrize: ", sum(becsum))

    Nspin = 1
    ia = 1 # atom index
    isp = Ham.atoms.atm2species[ia] # species index
    Nrmesh = Ham.pspots[isp].Nr
    l2 = (Ham.pspots[isp].lmax_rho + 1)^2

    # Calculate rho_lm
    rho_lm = zeros(Float64, Nrmesh, l2, Nspin)
    AE = true
    PAW_rho_lm!(AE, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm)
    println("sum rho_lm = ", sum(rho_lm))

    rho_rad = zeros(Float64, Nrmesh, Nspin)
    grad = zeros(Float64, Nrmesh, 3, Nspin) # gradient (r, ϕ, θ)
    grad2 = zeros(Float64, Nrmesh, Nspin) # square modulus of gradient

    atoms = Ham.atoms
    pspots = Ham.pspots
    pspotNL = Ham.pspotNL
    nx = pspotNL.paw.spheres[isp].nx

    if AE
        rho_core = pspots[isp].paw_data.ae_rho_atc
    else
        rho_core = pspots[isp].rho_atc
    end

    r = Ham.pspots[isp].r
    r2 = r.^2

    arho = zeros(Float64, Nrmesh) # 2nd dim removed 
    gradx = zeros(Float64, 3, Nrmesh) # 3rd dim removed
    grhoe2 = zeros(Float64, Nrmesh)


    # These arrays should depend on spin
    v1x = zeros(Float64, Nrmesh, Nspin)
    v2x = zeros(Float64, Nrmesh, Nspin)
    v1c = zeros(Float64, Nrmesh, Nspin)
    v2c = zeros(Float64, Nrmesh, Nspin) 
    # These arrays don't depend on spin
    sx = zeros(Float64, Nrmesh)
    sc = zeros(Float64, Nrmesh)

    e_rad = zeros(Float64, Nrmesh)
    gc_rad = zeros(Float64, Nrmesh, nx, Nspin)
    h_rad = zeros(Float64, Nrmesh, 3, nx, Nspin)

    energy = 0.0 # This should be accumulated for all ix
    spheres = pspotNL.paw.spheres

    #ix = 1
    #@assert ix <= nx

    for ix in 1:nx

        PAW_lm2rad!(ia, ix, atoms, pspots, pspotNL, rho_lm, rho_rad)
        PAW_gradient!(ia, ix, atoms, pspots, pspotNL,
            rho_lm, rho_rad, rho_core,
            grad2, grad
        )
        @printf("sum grad2 = %18.10e\n", sum(grad2))
        @printf("sum grad r     = %18.10e\n", sum(grad[:,1,:]))
        @printf("sum grad phi   = %18.10e\n", sum(grad[:,2,:]))
        @printf("sum grad theta = %18.10e\n", sum(grad[:,3,:]))

        for ir in 1:Nrmesh
            arho[ir,1] = abs(rho_rad[ir,1]/r2[ir] + rho_core[ir])
            @views gradx[1:3,ir] .= grad[ir,1:3] # transposed?
            grhoe2[ir] = grad[ir,1]^2 + grad[ir,2]^2 + grad[ir,3]^2
        end

        println("sum arho = ", sum(arho))
        println("sum gradx = ", sum(gradx))

        # sx, sc, v1x, v2x, v1c, v2c
        for ir in 1:Nrmesh
            sx[ir], v1x[ir,1], v2x[ir,1] = PWDFT.XC_x_pbe( arho[ir,1], grhoe2[ir] )
            sc[ir], v1c[ir,1], v2c[ir,1] = PWDFT.XC_c_pbe( arho[ir,1], grhoe2[ir] )
            #println("$ir $(arho[ir]) $(sx[ir])")
            # # This is needed in case of QE using Libxc
            # eex, vvx = PWDFT.XC_x_slater( arho[ir,1] )
            # eec, vvc = PWDFT.XC_c_pw( arho[ir,1] )
            # sx[ir] += eex*arho
            # sc[ir] += eec*arho
            # # Not yet for vvx and vvc
        end

        println("sum sx = ", sum(sx))
        println("sum sc = ", sum(sc))

        # radial stuffs
        for ir in 1:Nrmesh
            e_rad[ir] = (sx[ir] + sc[ir]) * r2[ir]
            gc_rad[ir,ix,1]  = ( v1x[ir,1] + v1c[ir,1] )
            @views h_rad[ir,1:3,ix,1] = ( v2x[ir,1] + v2c[ir,1] )*grad[ir,1:3,1]*r2[ir]
        end
    
        # integrate to obtain the energy
        energy += PWDFT.integ_simpson(Nrmesh, e_rad, pspots[isp].rab)*spheres[isp].ww[ix]

    end

    println("energy for all ix = ", energy)

end

main()