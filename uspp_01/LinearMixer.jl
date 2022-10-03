struct LinearMixer
    betamix::Float64
end

# Rhoe is not used
function LinearMixer(Rhoe::Matrix{Float64}, betamix)
    return LinearMixer(betamix)
end

# iterSCF is not used in this case
function do_mix!(mixer::LinearMixer, Rhoe, Rhoe_new, iterSCF)
    @views Rhoe[:] .= mixer.betamix*Rhoe_new[:] + (1 - mixer.betamix)*Rhoe[:]
    return
end

# Rhoe = β*Rhoe_new + (1-β)*Rhoe
#      = β*Rhoe_new + Rhoe - β*Rhoe
#      = β*(Rhoe_new - Rhoe) + Rhoe

function do_mix_precKerker!(mixer::LinearMixer, pw::PWGrid, Rhoe, Rhoe_new, iterSCF)
    dRhoe = precKerker(pw, Rhoe_new - Rhoe)
    @views Rhoe[:] .= mixer.betamix*dRhoe[:] + Rhoe[:]
    return
end

function precKerker(pw::PWGrid, R)
    Rg = R_to_G(pw, vec(R))
    idx_g2r = pw.gvec.idx_g2r
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    Ag = pw.gvec.G2_shells[2] # 1.0 # 
    for ig in 1:Ng
        ip = idx_g2r[ig]
        Rg[ip] = G2[ig]/(G2[ig] + Ag)*Rg[ip]
    end
    return real(G_to_R(pw, Rg))
end