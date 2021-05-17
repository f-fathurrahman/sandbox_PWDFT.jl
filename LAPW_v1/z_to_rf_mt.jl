function z_to_rf_mt!(
    mt_vars::MuffinTins,
    nr::Int64, nri::Int64,
    zfmt, rfmt
)
    #
    lmmaxi = mt_vars.lmmaxi
    lmaxi = mt_vars.lmaxi
    #
    lmaxo = mt_vars.lmaxo
    lmmaxo = mt_vars.lmmaxo
    #
    i = 1
    for ir in 1:nri
        @views z_to_rf_lm!(lmaxi, zfmt[i:end], rfmt[i:end])
        i = i + lmmaxi
    end
    for ir in (nri+1):nr
        @views z_to_rf_lm!(lmaxo, zfmt[i:end], rfmt[i:end])
        i = i + lmmaxo
    end
    return
end