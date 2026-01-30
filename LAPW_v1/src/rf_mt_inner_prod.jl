
# originally rfmtinp
function rf_mt_inner_prod(isp, mt_vars, rfmt1, rfmt2)

    lmaxi = mt_vars.lmaxi
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    nr = mt_vars.nrmt[isp]
    nri = mt_vars.nrmti[isp]
    wr = mt_vars.wrmt[isp]

    fr = zeros(Float64, nr)
    i = 1
    # inner part of muffin-tin
    if lmaxi == 1
        for ir in  1:nri
            fr[ir] = rfmt1[i]*rfmt2[i] + rfmt1[i+1]*rfmt2[i+1] +
                     rfmt1[i+2]*rfmt2[i+2] + rfmt1[i+3]*rfmt2[i+3]
            i += 4 # increment
        end
    else
        n = lmmaxi - 1
        for ir in 1:nri
            fr[ir] = dot(rfmt1[i:i+n], rfmt2[i:i+n])
            i += lmmaxi # increment
        end
    end
    # outer part of muffin-tin
    n = lmmaxo - 1
    for ir in (nri+1):nr
        fr[ir] = dot( rfmt1[i:i+n], rfmt2[i:i+n])
        i += lmmaxo
    end
    # integrate
    return dot(wr, fr)
end
