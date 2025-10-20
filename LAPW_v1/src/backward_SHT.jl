# From rbsht
function backward_SHT!(
    mt_vars, isp,
    rfmt_in::Vector{Float64},
    rfmt_out::Vector{Float64}
)
    nr = mt_vars.nrmt[isp]
    nri = mt_vars.nrmti[isp]
    nro = nr - nri

    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo

    idx_inner = 1:lmmaxi*nri
    RB = mt_vars.SHT.rbshti
    @views V_in = reshape(rfmt_in[idx_inner], (lmmaxi,nri))
    @views V_out = reshape(rfmt_out[idx_inner], (lmmaxi,nri))
    @views V_out[:,:] = RB[:,:] * V_in[:,:]

    idx_outer = (lmmaxi*nri + 1):(lmmaxi*nri + lmmaxo*nro)
    RB = mt_vars.SHT.rbshto
    @views V_in = reshape(rfmt_in[idx_outer], (lmmaxo,nro))
    @views V_out = reshape(rfmt_out[idx_outer], (lmmaxo,nro))
    @views V_out[:,:] = RB[:,:] * V_in[:,:]

    return
end

# From zbsht (ComplexF64 version)
function backward_SHT!(
    mt_vars, isp,
    zfmt_in::Vector{ComplexF64},
    zfmt_out::Vector{ComplexF64}
)
    nr = mt_vars.nrmt[isp]
    nri = mt_vars.nrmti[isp]
    nro = nr - nri

    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo

    idx_inner = 1:lmmaxi*nri
    ZB = mt_vars.SHT.zbshti
    @views V_in = reshape(zfmt_in[idx_inner], (lmmaxi,nri))
    @views V_out = reshape(zfmt_out[idx_inner], (lmmaxi,nri))
    @views V_out[:,:] = ZB[:,:] * V_in[:,:]

    idx_outer = (lmmaxi*nri + 1):(lmmaxi*nri + lmmaxo*nro)
    ZB = mt_vars.SHT.zbshto
    @views V_in = reshape(zfmt_in[idx_outer], (lmmaxo,nro))
    @views V_out = reshape(zfmt_out[idx_outer], (lmmaxo,nro))
    @views V_out[:,:] = ZB[:,:] * V_in[:,:]

    return
end

