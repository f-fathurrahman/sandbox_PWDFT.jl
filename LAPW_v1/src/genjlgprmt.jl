function genjlgprmt!(rmt, lmax, G2, jlgprmt)
    Ng = size(G2,1)
    Nspecies = size(rmt,1)
    for isp in 1:Nspecies
        r = rmt[isp]
        for ig in 1:Ng
            t1 = sqrt(G2[ig])*r
            @views sbessel!( lmax, t1, jlgprmt[:,ig,isp] )
        end
    end
    return
end
