mutable struct AndersonMixer
    betamix::Float64
    mixdim::Int64
    df::Matrix{Float64}
    dv::Matrix{Float64}
end


function AndersonMixer(Rhoe::Matrix{Float64}, betamix; mixdim=8)
    df = zeros(Float64, length(Rhoe), mixdim)
    dv = zeros(Float64, length(Rhoe), mixdim)
    return AndersonMixer(betamix, mixdim, df, dv)
end


function do_mix!(
    mixer::AndersonMixer,
    Rhoe, Rhoe_new,
    iterSCF
)
    mix_anderson!( Rhoe, Rhoe_new,
        mixer.betamix, mixer.df, mixer.dv,
        iterSCF, mixer.mixdim
    )
    return
end

#
# This function is adapted from Anderson mixing function in KSSOLV
#
# vin is the result
# vin will be overwritten for next SCF iteration
function _do_mix_anderson!(
    vin, vout,
    betamix::Float64, df::Array{Float64,2}, dv::Array{Float64,2},
    iter::Int64, mixdim::Int64
)
    # Residual
    dvout = vout[:] - vin[:]

    iterused = min(iter-1,mixdim)
    ipos = iter - 1 - floor(Int64, (iter-2)/mixdim)*mixdim
    println("mix_anderson: ipos = ", ipos)

    if iter > 1
        @views df[:,ipos] = df[:,ipos] - dvout[:]
        @views dv[:,ipos] = dv[:,ipos] - vin[:]
    end

    vinsave  = copy(vin)
    dvoutsave = copy(dvout)

    if iter > 1
        gammas = pinv(df[:,1:iterused])*dvout  
        for i in 1:iterused
            @views vin[:]  = vin[:] - gammas[i] * dv[:,i]
            @views dvout[:] = dvout[:] - gammas[i] * df[:,i]
        end
    end

    inext = iter - floor( Int64, (iter - 1) / mixdim) * mixdim

    @views df[:,inext] = dvoutsave[:]
    @views dv[:,inext] = vinsave[:]

    @views vin[:] = vin[:] + betamix*dvout[:]

    return
end

