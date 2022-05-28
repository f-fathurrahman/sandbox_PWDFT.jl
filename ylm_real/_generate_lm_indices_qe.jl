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