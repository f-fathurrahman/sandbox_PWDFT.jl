function rotrfmt!(mt_vars, isp::Int64, R, rfmt1, rfmt2; coarse=false)
    # XXX It seems that `isp` can be mistaken as `ia` (?)
    if coarse
        nr = mt_vars.nrcmt[isp]
        nri = mt_vars.nrcmti[isp]        
    else
        nr = mt_vars.nrmt[isp]
        nri = mt_vars.nrmti[isp]
    end
    nro = nr - nri

    # inner part of muffin-tin
    lmaxi = mt_vars.lmaxi
    lmmaxi = mt_vars.lmmaxi
    idx_inner = 1:lmmaxi*nri
    @views rotrflm!(R, lmaxi, nri, lmmaxi, rfmt1[idx_inner], rfmt2[idx_inner])
    
    # outer part of muffin-tin
    lmaxo = mt_vars.lmaxo
    lmmaxo = mt_vars.lmmaxo
    idx_outer = (lmmaxi*nri + 1):(lmmaxi*nri + lmmaxo*nro)
    @views rotrflm!(R, lmaxo, nro, lmmaxo, rfmt1[idx_outer], rfmt2[idx_outer])
    
    return
end
