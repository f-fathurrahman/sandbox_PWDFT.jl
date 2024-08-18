function hmlistl!(ik, pw, cfunig, vsig, nmat, H)
    #
    Ngwk = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g
    idx_g2miller = pw.gvec.idx_g2miller
    G = pw.gvec.G
    @views kvec = pw.gvecw.kpoints.k[:,ik]
    Gk_j = zeros(Float64, 3)
    Gk_i = zeros(Float64, 3) 
    #
    for j in 1:Ngwk
        k = (j-1)*nmat[ik]
        jg = idx_gw2g[ik][j]
        Gk_j[:] .= G[:,jg] .+ kvec[:]
        j_idx = idx_g2miller[jg] # tuple
        for i in 1:j
            k += 1
            ig = idx_gw2g[ik][i]
            Gk_i[:] .= G[:,ig] .+ kvec[:]
            t1 = 0.5*dot(Gk_i, Gk_j)
            ij_idx = idx_g2miller[ig] .- j_idx # tuple
            #
            ijg = findfirst( isequal(ij_idx), idx_g2miller )
            H[k] += t1*cfunig[ijg] + vsig[ijg] 
        end 
    end
    return
end

