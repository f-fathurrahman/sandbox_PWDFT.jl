using Printf
using OffsetArrays

include("Ylm_real_qe.jl")

# Using QE's convention
function _generate_lm_indices_qe(lmax)
    lmmax = (lmax+1)^2
    # index to (l,m) pairs
    idxlm = OffsetArray( zeros(Int64,lmax+1, 2*lmax+1),
        0:lmax,-lmax:lmax)
    idxil = zeros(Int64, lmmax)
    idxim = zeros(Int64, lmmax)
    lm = 0
    for l in 0:lmax
        # m = 0
        lm = lm + 1
        idxlm[l,0] = lm
        idxil[lm] = l
        idxim[lm] = 0
        for m in 1:l
            lm = lm + 1
            idxlm[l,m] = lm
            idxil[lm] = l
            idxim[lm] = m
            # negative m
            lm = lm + 1
            idxlm[l,-m] = lm
            idxil[lm] = l
            idxim[lm] = -m

        end
    end
    return idxlm, idxil, idxim
end

function test_main()

    g = [0.1, 1.0, 2.0]
    Ng = 1
    lmax = 3
    lmmax = (lmax+1)^2

    idxlm, idxil, idxim = _generate_lm_indices_qe(lmax)
    
    Ylm = zeros(lmmax)
    Ylm_real_qe!(lmax, g, Ylm)
    
    #for l in 0:lmax
    #    println()
    #    for m in -l:l
    #        lm = idxlm[l,m]
    #        @printf("l=%3d m=%3d lm=%3d Ylm = %18.10f\n", l, m, lm, Ylm[lm])
    #    end
    #end

    println()
    for lm in 1:lmmax
        l = idxil[lm]
        m = idxim[lm]
        @printf("l=%3d m=%3d lm=%3d Ylm = %18.10f\n", l, m, lm, Ylm[lm])
    end

end

test_main()
