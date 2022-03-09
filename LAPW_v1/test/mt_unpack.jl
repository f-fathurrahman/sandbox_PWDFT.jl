function mt_unpack!(
  lmmaxo::Int64, lmmaxi::Int64,
  nr::Int64, nri::Int64,
  rfmt1, rfmt2
)
    i = 1
    j = 1
    for ir in 1:nri
        @views rfmt2[j:j+lmmaxi-1] = rfmt1[i:i+lmmaxi-1]
        i = i + lmmaxi
        k = j + lmmaxi
        j = j + lmmaxo
        @views rfmt2[k:j-1] .= 0.0
    end
    k = lmmaxo*(nr - nri)
    rfmt2[j:j+k-1] = rfmt1[i:i+k-1]
    return 
end

