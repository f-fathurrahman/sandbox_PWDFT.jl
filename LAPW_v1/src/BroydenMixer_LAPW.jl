mutable struct BroydenMixer_LAPW
    betamix::Float64
    mixdim::Int64
    df_ir::Matrix{Float64}
    dv_ir::Matrix{Float64}
    df_mt::Vector{Matrix{Float64}}
    dv_mt::Vector{Matrix{Float64}}
end


function BroydenMixer_LAPW(
    rho_mt::Union{Vector{Vector{Float64}}, Vector{Matrix{Float64}}},
    rho_ir::Union{Vector{Float64}, Matrix{Float64}},
    betamix;
    mixdim=8
)
    df_ir = zeros(Float64, length(rho_ir), mixdim)
    dv_ir = zeros(Float64, length(rho_ir), mixdim)

    Natoms = size(rho_mt, 1)
    df_mt = Vector{Matrix{Float64}}(undef, Natoms)
    dv_mt = Vector{Matrix{Float64}}(undef, Natoms)
    for ia in 1:Natoms
        df_mt[ia] = zeros(Float64, length(rho_mt[ia]), mixdim)
        dv_mt[ia] = zeros(Float64, length(rho_mt[ia]), mixdim)
    end
    return BroydenMixer_LAPW(betamix, mixdim, df_ir, dv_ir, df_mt, dv_mt)
end


function do_mix_LAPW!(
    mixer::BroydenMixer_LAPW,
    rho_ir_in, rho_ir_out_,
    rho_mt_in, rho_mt_out_,
    iterSCF::Int64
)
    # XXX: alphamix -> mixer.betamix for now
    PWDFT._do_mix_broyden!(
        rho_ir_in, rho_ir_out_,
        mixer.betamix,
        iterSCF, mixer.mixdim,
        mixer.df_ir, mixer.dv_ir
    )
    Natoms = size(rho_mt_in, 1)
    for ia in 1:Natoms
        PWDFT._do_mix_broyden!(
            rho_mt_in[ia], rho_mt_out_[ia],
            mixer.betamix,
            iterSCF, mixer.mixdim,
            mixer.df_mt[ia], mixer.dv_mt[ia]
        )
    end
    return
end
