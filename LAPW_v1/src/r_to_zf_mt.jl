function r_to_zf_mt!(
    mt_vars,
    nr::Int64, nri::Int64,
    rfmt, zfmt
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
        @views r_to_zf_lm!( lmaxi, rfmt[i:end], zfmt[i:end] )
        i = i + lmmaxi
    end
    for ir in nri+1:nr
        @views r_to_zf_lm( lmaxo, rfmt[i:end], zfmt[i:end] )
        i = i + lmmaxo
    end
    return
end

