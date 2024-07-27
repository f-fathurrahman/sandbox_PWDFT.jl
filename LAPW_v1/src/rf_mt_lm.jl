function rf_mt_lm!(lm, isp, mt_vars, rfmt, fr)

    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    nri = mt_vars.nrmti[isp]
    nr = mt_vars.nrmt[isp]

    if lm > lmmaxi
        fr[1:nri] .= 0.0
    else 
        i = lm
        for ir in 1:nri
            fr[ir] = rfmt[i]
            i = i + lmmaxi
        end
    end
  
    iro = nri + 1
    if lm > lmmaxo
        fr[iro:nr] .= 0.0
    else
        i = lmmaxi*nri + lm
        for ir in iro:nr
            fr[ir] = rfmt[i]
            i = i + lmmaxo
        end
    end

    return
end
