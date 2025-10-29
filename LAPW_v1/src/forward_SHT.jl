# From rfsht
function forward_SHT!(
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
    RF = mt_vars.SHT.rfshti
    @views V_in = reshape(rfmt_in[idx_inner], (lmmaxi,nri))
    @views V_out = reshape(rfmt_out[idx_inner], (lmmaxi,nri))
    @views V_out[:,:] = RF[:,:] * V_in[:,:]

    idx_outer = (lmmaxi*nri + 1):(lmmaxi*nri + lmmaxo*nro)
    RF = mt_vars.SHT.rfshto
    @views V_in = reshape(rfmt_in[idx_outer], (lmmaxo,nro))
    @views V_out = reshape(rfmt_out[idx_outer], (lmmaxo,nro))
    @views V_out[:,:] = RF[:,:] * V_in[:,:]

    return
end

# From zfsht (ComplexF64 version)
function forward_SHT!(
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
    ZF = mt_vars.SHT.zfshti
    @views V_in = reshape(zfmt_in[idx_inner], (lmmaxi,nri))
    @views V_out = reshape(zfmt_out[idx_inner], (lmmaxi,nri))
    @views V_out[:,:] = ZF[:,:] * V_in[:,:]

    idx_outer = (lmmaxi*nri + 1):(lmmaxi*nri + lmmaxo*nro)
    ZF = mt_vars.SHT.zfshto
    @views V_in = reshape(zfmt_in[idx_outer], (lmmaxo,nro))
    @views V_out = reshape(zfmt_out[idx_outer], (lmmaxo,nro))
    @views V_out[:,:] = ZF[:,:] * V_in[:,:]

    return
end

