function genffacgp(pw, rmt; gvec_full=nothing)
    if isnothing(gvec_full)
        gvec = pw.gvec
    else
        gvec = gvec_full
    end
    Nspecies = size(rmt, 1)
    ffacgp = zeros(Float64, gvec.Ng, Nspecies)
    genffacgp!(pw, rmt, ffacgp, gvec_full=gvec_full)
    return ffacgp
end

function genffacgp!(pw::PWGrid, rmt::Vector{Float64}, ffacgp; gvec_full=nothing)
    if isnothing(gvec_full)
        gvec = pw.gvec
    else
        gvec = gvec_full
    end

    G2 = gvec.G2
    Ng = gvec.Ng
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
