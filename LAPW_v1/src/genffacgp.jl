function genffacgp!(pw::PWGrid, rmt::Vector{Float64}, ffacgp)
    # Only G2 is considered, so gpc => pw.gvec.G2

    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    CellVolume = pw.CellVolume
    t1 = 4π / CellVolume
    
    Nspecies = size(rmt, 1)
    for isp in 1:Nspecies
        # First G-vector is 0
        ffacgp[1,isp] = (t1/3.0)*rmt[isp]^3
        # XXX Need to check from small G2[ig] ?
        for ig in 2:Ng
            G2_len = sqrt(G2[ig])
            t2 = G2_len*rmt[isp]
            ffacgp[ig,isp] = t1*( sin(t2) - t2*cos(t2) )/( G2_len^3 )
        end
    end 
    return
end
