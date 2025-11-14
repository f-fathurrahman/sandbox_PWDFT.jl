function zfmtinp(mt_vars, isp, wr, zfmt1, zfmt2; coarse=false)

    if coarse
        nr = mt_vars.nrcmt[isp]
        nri = mt_vars.nrcmti[isp]
    else
        nr = mt_vars.nrmt[isp]
        nri = mt_vars.nrmti[isp]
    end

    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    lmaxi = mt_vars.lmaxi

    fr1 = zeros(Float64, nr)
    fr2 = zeros(Float64, nr)
    # compute the dot-products for each radial point
    i = 1
    if lmaxi == 1
        for ir in 1:nri
            z1 = conj(zfmt1[i])*zfmt2[i] + conj(zfmt1[i+1])*zfmt2[i+1] +
                 conj(zfmt1[i+2])*zfmt2[i+2] + conj(zfmt1[i+3])*zfmt2[i+3]
            fr1[ir] = real(z1)
            fr2[ir] = imag(z1)
            i += 4  # HARCODED?? lradstp
        end
    else
        for ir in 1:nri
            z1 = BLAS.dotc(lmmaxi, zfmt1[i:i+lmmaxi-1], 1 , zfmt2[i:i+lmmaxi-1], 1)
            fr1[ir] = real(z1)
            fr2[ir] = imag(z1)
            i += lmmaxi
        end
    end
    #
    for ir in (nri+1):nr
        z1 = BLAS.dotc(lmmaxo, zfmt1[i:i+lmmaxo-1], 1, zfmt2[i:i+lmmaxo-1], 1)
        fr1[ir] = real(z1)
        fr2[ir] = imag(z1)
        i += lmmaxo
    end
    # integrate over r
    return dot(wr, fr1) + im*dot(wr, fr2)
end