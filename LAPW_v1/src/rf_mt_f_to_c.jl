function rf_mt_f_to_c!(
    mt_vars,
    nrc ,nrci, rfmt, rfcmt
)
    lmmaxo = mt_vars.lmmaxo
    lmmaxi = mt_vars.lmmaxi
    lradstp = mt_vars.lradstp
    i = 1
    j = 1
    n = lmmaxi*lradstp
    for irc in 1:nrci
        rfcmt[j:(j+lmmaxi-1)] = rfmt[i:(i+lmmaxi-1)]
        i = i + n
        j = j + lmmaxi
    end
    i = i + (lradstp-1)*(lmmaxo-lmmaxi)
    n = lmmaxo*lradstp
    for irc in (nrci+1):nrc
        rfcmt[j:(j+lmmaxo-1)] = rfmt[i:(i+lmmaxo-1)]
        i = i + n
        j = j + lmmaxo
    end
    return
end