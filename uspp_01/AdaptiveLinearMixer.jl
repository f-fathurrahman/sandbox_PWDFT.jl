mutable struct AdaptiveLinearMixer
    betamix::Float64
    MAXBETA::Float64
    betav::Vector{Float64}
    f::Vector{Float64}
end

function AdaptiveLinearMixer(Rhoe, betamix; MAXBETA=0.8)
    betav = betamix*ones(Float64, length(Rhoe))
    f = zeros(Float64, length(Rhoe))
    return AdaptiveLinearMixer(betamix, MAXBETA, betav, f)
end


function do_mix!(
    mixer::AdaptiveLinearMixer,
    Rhoe, Rhoe_new,
    iterSCF
)
    _do_mix_adaptive!(Rhoe, Rhoe_new, mixer.betamix, mixer.betav, mixer.f, MAXBETA=mixer.MAXBETA)
    return
end



function _do_mix_adaptive!( mu, nu, betamix::Float64, betav, f; MAXBETA=0.8 )
    Npts = length(mu)
    for i in 1:Npts
        t = nu[i] - mu[i]
        if t*f[i] >= 0.0
            betav[i] = betav[i] + betamix
            if betav[i] > MAXBETA
                betav[i] = MAXBETA
            end
        else
            betav[i] = 0.5*( betav[i] + betamix )
        end        
        f[i] = t
        mu[i] = betav[i]*nu[i] + ( 1.0 - betav[i] )*mu[i]
    end
    return
end