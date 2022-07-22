function dense_to_smooth!(
    pw::PWGrid, v_in, v_out
)

    # Immediate return if not using dual grid
    if !pw.using_dual_grid
        v_out[:] = v_in[:]
        return
    end

    NptsDense = prod(pw.Ns)
    NptsSmooth = prod(pw.Nss) # not used?

    aux_in = zeros(ComplexF64, NptsDense)

    # Copy all data to aux_in
    aux_in[1:NptsDense] .= v_in[1:NptsDense]

    R_to_G!(pw, aux_in) # need a way to tell which plan to use

    #println("Some aux_in after R_to_G! (scaled)")
    #for i in 1:10
    #    @printf("%5d %18.10f %18.10f\n", i, aux_in[i].re/NptsDense, aux_in[i].im/NptsDense)
    #end

    #println("sum aux_in, scaled = ", sum(aux_in)/NptsDense)

    # Zero out aux_out (not needed)

    aux_out = zeros(ComplexF64, NptsSmooth)
    Ngs = pw.gvecs.Ng
    #println("Ngs = ", Ngs)
    for ig in 1:Ngs
        ip_out = pw.gvecs.idx_g2r[ig] # smooth
        ip_in  = pw.gvec.idx_g2r[ig] # dense
        aux_out[ip_out] = aux_in[ip_in]
    end

    #println("Compare ip_in, ip_out")
    #for ig in 1:10
    #    ip_out = pw.gvecs.idx_g2r[ig]
    #    ip_in  = pw.gvec.idx_g2r[ig]
    #    @printf("%5d %5d\n", ip_in, ip_out)
    #end

    #println("Some aux_out before G_to_R! (scaled)")
    #for i in NptsSmooth-10:NptsSmooth
    #    @printf("%8d %18.10f %18.10f\n", i, aux_out[i].re/NptsDense, aux_out[i].im/NptsDense)
    #end

    #println("Before G_to_R!: sum(aux_out), scaled = ", sum(aux_out)/NptsDense)
    
    G_to_R!(pw, aux_out, smooth=true)
    # FIXME: Manual normalization
    aux_out *= (NptsSmooth/NptsDense)  # XXX This is important!

    #println("Some aux_out after G_to_R!")
    #for i in 1:10
    #    @printf("%5d %18.10f %18.10f\n", i, aux_out[i].re, aux_out[i].im)
    #end

    #println("After G_to_R!: sum(aux_out) = ", sum(aux_out))

    v_out[1:NptsSmooth] = real(aux_out[1:NptsSmooth])

    return
end