function genffacgp!(pw::PWGrid, isp::Int64, rmt_isp::Float64, ffacgp)
    # Only G2 is considered, so gpc => pw.gvec.G2

    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    CellVolume = pw.CellVolume
    t1 = 4π / CellVolume
    
    # First G-vector is 0
    t2 = G2[1]*rmt_isp
    ffacgp[1] = t1*( sin(t2) - t2*cos(t2) )/( G2[ig]^3 )
    # In the original code this is done for G2 < epslat

    for ig in 2:Ng
        ffacgp[ig] = (t1/3.0)*rmt_isp^3
    end 
    return
end
