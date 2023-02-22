using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("PAWAtomicSphere.jl")
#include("dYlm_real_qe.jl")

function main()
    Ham, pwinput = init_Ham_from_pwinput()

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    # Some constants from paw_variables.f90
    lm_fact = 3   # To converge E_xc integrate up to LM = lm_fact * lm_max
    lm_fact_x = 3 # As above, for gradient corrected functionals
    xlm = 2       # Additional factor to add to have a good grad.corr.
    radial_grad_style = 0 # = 0 or 1, algorithm to use for d/dr

    is_gga = false
    lmax_safe = 0
    lmax_add = 0

    paw_spheres = Vector{PAWAtomicSphere}(undef,Nspecies)

    for isp in 1:Nspecies
        
        psp = Ham.pspots[isp]

        if !psp.is_paw
            continue
        end
        
        if psp.lmax_rho == 0
            # no need for more than one direction, when it is spherical!
            lmax_safe = 0
            lmax_add  = 0
        else
            if is_gga
                # Integrate up to a higher maximum lm if using gradient
                # correction check expression for d(y_lm)/d\theta for details
                lmax_safe = lm_fact_x * psp.lmax_rho
                lmax_add  = xlm
            else
                # no gradient correction:
                lmax_safe = lm_fact * psp.lmax_rho
                lmax_add  = 0 
            end
        end

        println("lmax_safe = ", lmax_safe)
        println("lmax_add  = ", lmax_add)

        paw_spheres[isp] = PAWAtomicSphere(lmax_safe, lmax_add)

    end


end

main()