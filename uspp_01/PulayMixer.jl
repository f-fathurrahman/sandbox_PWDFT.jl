mutable struct PulayMixer
    betamix::Float64
    mixdim::Int64
    XX::Matrix{Float64}
    FF::Matrix{Float64}
    x_old::Matrix{Float64}
    f_old::Matrix{Float64}
end


function PulayMixer(Rhoe, betamix; mixdim=4)
    XX = zeros(Float64, length(Rhoe), mixdim)
    FF = zeros(Float64, length(Rhoe), mixdim)
    x_old = zeros(Float64, size(Rhoe))
    f_old = zeros(Float64, size(Rhoe))
    return PulayMixer(betamix, mixdim, XX, FF, x_old, f_old)
end


function do_mix!(
    mixer::PulayMixer,
    Rhoe, Rhoe_new,
    iterSCF
)
    _do_mix_pulay!( Rhoe, Rhoe_new, mixer.betamix,
        mixer.XX, mixer.FF,
        iterSCF, mixer.mixdim,
        mixer.x_old, mixer.f_old
    )
    return
end


function _do_mix_pulay!(
    x, gx, beta, X, F, iter, MIXDIM::Int64,
    x_old, f_old
)

    f = gx[:] - x[:]

    if iter == 1
        dx = copy(x)
        df = copy(f)
    else
        dx = x[:] - x_old[:]
        df = f[:] - f_old[:]
    end

    # shift history
    if iter <= MIXDIM
        X[:,iter] = dx[:]
        F[:,iter] = df[:]
    else
        for i in 1:MIXDIM-1
            X[:,i] = X[:,i+1]
            F[:,i] = F[:,i+1]
        end
        # new history, put at the end
        X[:,MIXDIM] = dx[:]
        F[:,MIXDIM] = df[:]
    end

    x_old[:] = x[:]
    f_old[:] = f[:]

    # Pulay mixing begins at iter MIXDIM+1
    if iter > MIXDIM
        addv = (X + beta*F)*inv(F'*F)*(F'*f)
        x[:] = x[:] + beta*f - addv
    else
        x[:] = x[:] + beta*f
    end
    
    return

end


function do_mix_precKerker!(
    mixer::PulayMixer,
    pw::PWGrid,
    Rhoe, Rhoe_new,
    iterSCF
)
    _do_mix_pulay!( pw, Rhoe, Rhoe_new, mixer.betamix,
        mixer.XX, mixer.FF,
        iterSCF, mixer.mixdim,
        mixer.x_old, mixer.f_old
    )
    return
end


function _do_mix_pulay!(
    pw::PWGrid,
    x, gx, beta, X, F, iter, MIXDIM::Int64,
    x_old, f_old
)
    f = precKerker(pw, gx[:] - x[:])
    if iter == 1
        dx = copy(x)
        df = copy(f)
    else
        dx = x[:] - x_old[:]
        df = f[:] - f_old[:]
    end
    # shift history
    if iter <= MIXDIM
        X[:,iter] = dx[:]
        F[:,iter] = df[:]
    else
        for i in 1:MIXDIM-1
            X[:,i] = X[:,i+1]
            F[:,i] = F[:,i+1]
        end
        # new history, put at the end
        X[:,MIXDIM] = dx[:]
        F[:,MIXDIM] = df[:]
    end
    x_old[:] = x[:]
    f_old[:] = f[:]
    # Pulay mixing begins at iter MIXDIM+1
    if iter > MIXDIM
        addv = (X + beta*F)*inv(F'*F)*(F'*f)
        x[:] = x[:] + beta*f - addv
    else
        x[:] = x[:] + beta*f
    end
    return
end
