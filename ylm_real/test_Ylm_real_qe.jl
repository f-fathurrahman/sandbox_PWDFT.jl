using Printf
using OffsetArrays
using PWDFT: Ylm_real

include("Ylm_real_qe.jl")
include("_generate_lm_indices_qe.jl")

function test_main()

    g = [0.1, 1.0, 2.0]
    Ng = 1
    lmax = 3
    lmmax = (lmax+1)^2

    idxlm, idxil, idxim = _generate_lm_indices_qe(lmax)
    
    Ylm = zeros(lmmax)
    Ylm_real_qe!(lmax, g, Ylm)
    
    # Using the "natural" loop (sum over l (sum over m))
    #for l in 0:lmax
    #    println()
    #    for m in -l:l
    #        lm = idxlm[l,m]
    #        @printf("l=%3d m=%3d lm=%3d Ylm = %18.10f\n", l, m, lm, Ylm[lm])
    #    end
    #end

    # Using linear indexing (combined l,m index), starting from 1 
    println()
    for lm in 1:lmmax
        l = idxil[lm]
        m = idxim[lm]
        yy = Ylm_real(l, m, g) # * (-1)^m (Need this factor to be the same as QE?)
        Δ = Ylm[lm] - yy
        @printf("l=%3d m=%3d lm=%3d Ylm = %18.10f %18.10f %18.10f\n",
            l, m, lm, Ylm[lm], yy, Δ)
    end

end

test_main()
