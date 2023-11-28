using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))


function my_PAW_rad2lm!(
    ia,
    atoms, pspotNL,
    lmax_loc::Int64,
    F_rad, F_lm
)
    Nspin = size(F_rad, 3)
    Nrmesh = size(F_rad, 1)
    isp = atoms.atm2species[ia]
    sphere = pspotNL.paw.spheres[isp]

    for ispin in 1:Nspin, lm in 1:lmax_loc^2
        @views F_lm[:,lm,ispin] .= 0.0
        for ix in 1:sphere.nx, j in 1:Nrmesh
            F_lm[j,lm,ispin] = F_lm[j,lm,ispin] + F_rad[j,ix,ispin] * sphere.wwylm[ix,lm]
        end
    end

    file = open("fort.222julia", "w")
    for ispin in 1:Nspin
        for lm in 1:lmax_loc^2
            for j in 1:Nrmesh
                @printf(" %5d%5d%18.10f\n", j, lm, F_lm[j,lm,ispin])
            end
        end
    end

    return
end

function main(;filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    ia = 1 # choose atom index

    Nspin = 1
    atoms = Ham.atoms
    pspotNL = Ham.pspotNL
    isp = atoms.atm2species[ia]
    Nrmesh = Ham.pspots[isp].Nr
    l2 = (Ham.pspots[isp].lmax_rho + 1)^2
    spheres = Ham.pspotNL.paw.spheres
    nx = spheres[isp].nx
    lmax_loc = Ham.pspots[isp].lmax_rho + 1

    gc_rad = zeros(Float64, Nrmesh, nx, Nspin)
    gc_lm = zeros(Float64, Nrmesh, l2, Nspin)

    gc_rad[:,:,1] .= 1.1

    println("Nrmesh = ", Nrmesh)

    PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc, gc_rad, gc_lm)

    println("lmax_loc = ", lmax_loc)
    println("sum wwylm = ", sum(spheres[isp].wwylm))
    println("sum gc_rad = ", sum(gc_rad))
    println("sum gc_lm = ", sum(gc_lm))
end

main()
