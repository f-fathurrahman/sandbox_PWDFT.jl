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
    AE = false
    PAW_rho_lm!(AE, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm)
    println("sum rho_lm = ", sum(rho_lm))

    rho_rad = zeros(Float64, Nrmesh, Nspin)
    grad = zeros(Float64, Nrmesh,3, Nspin) # gradient (r, ϕ, θ)
    grad2 = zeros(Float64, Nrmesh, Nspin) # square modulus of gradient

    atoms = Ham.atoms
    pspots = Ham.pspots
    pspotNL = Ham.pspotNL

    if AE
        rho_core = pspots[isp].paw_data.ae_rho_atc
    else
        rho_core = pspots[isp].rho_atc
    end

    ix = 1
    PAW_lm2rad!(ia, ix, atoms, pspots, pspotNL, rho_lm, rho_rad)
    PAW_gradient!(ia, ix, atoms, pspots, pspotNL,
        rho_lm, rho_rad, rho_core,
        grad2, grad
    )
    @printf("sum grad2 = %18.10e\n", sum(grad2))
    @printf("sum grad r     = %18.10e\n", sum(grad[:,1,:]))
    @printf("sum grad phi   = %18.10e\n", sum(grad[:,2,:]))
    @printf("sum grad theta = %18.10e\n", sum(grad[:,3,:]))
end

main()