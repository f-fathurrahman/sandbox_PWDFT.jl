struct LinearMixer
    betamix::Float64
end

# Rhoe is not used
function LinearMixer(Rhoe::Matrix{Float64}, betamix)
    return LinearMixer(betamix)
end

# iterSCF is not used in this case
function do_mix!(mixer, Rhoe, Rhoe_new, iterSCF)
    @views Rhoe[:] .= mixer.betamix*Rhoe_new[:] + (1 - mixer.betamix)*Rhoe[:]
    return
end