function olpistl!(ik, pw, cfunig, O)
    nmatk = size(O, 1)
    Ngwk = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g
    idx_g2miller = pw.gvec.idx_g2miller
    #
    for j in 1:Ngwk
        k = (j-1)*nmatk
        jg = idx_gw2g[ik][j]
        j_idx = idx_g2miller[jg] # tuple
        for i in 1:j
            k += 1
            ig = idx_gw2g[ik][i]
            ij_idx = idx_g2miller[ig] .- j_idx # tuple
            ijg = findfirst( isequal(ij_idx), idx_g2miller )
            O[k] += cfunig[ijg]
        end 
    end
    return
end

