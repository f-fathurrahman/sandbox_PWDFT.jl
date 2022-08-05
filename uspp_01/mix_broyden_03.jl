function mix_broyden_03!(
    Ham::Hamiltonian, 
    Rhoe_in, Rhoe_out,
    alphamix::Float64, iter::Int64, n_iter::Int64,
    df, dv
)
    pw = Ham.pw
    Nelectrons = Ham.electrons.Nelectrons

    ctmp_in = zeros(ComplexF64, length(Rhoe_in))
    ctmp_out = zeros(ComplexF64, length(Rhoe_out))

    ctmp_in[:] = Rhoe_in[:]
    ctmp_out[:] = Rhoe_out[:]

    R_to_G!(pw, ctmp_in)
    R_to_G!(pw, ctmp_out)

    # FIXME: only works for Nspin = 1
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    RhoeG_in = zeros(ComplexF64, Ng)
    RhoeG_out = zeros(ComplexF64, Ng)

    for ig in 1:Ng
        ip = idx_g2r[ig]
        RhoeG_in[ig] = ctmp_in[ip]
        RhoeG_out[ig] = ctmp_out[ip]
    end

    mix_broyden_03!(pw, Nelectrons, RhoeG_in, RhoeG_out, alphamix, iter, n_iter, df, dv)

    for ig in 1:Ng
        ip = idx_g2r[ig]
        ctmp_in[ip] = RhoeG_in[ig]
        ctmp_out[ip] = RhoeG_out[ig]
    end

    # Back to R-space
    G_to_R!(pw, ctmp_in)
    G_to_R!(pw, ctmp_out)

    Rhoe_in[:] = real(ctmp_in[:])
    Rhoe_out[:] = real(ctmp_out[:])

    return
end



#
# Adapted from PWSCF/Yambo
#
function mix_broyden_03!( pw, Nelectrons,
    deltain, deltaout_, alphamix::Float64,
    iter::Int64, n_iter::Int64,
    df, dv
)
    
    deltaout = copy(deltaout_)  # do not replace deltaout_

    maxter = 8
    wg0 = 0.01
    wg = ones(maxter)

    deltainsave = copy( deltain )
    #
    # iter_used = iter-1  IF iter <= n_iter
    # iter_used = n_iter  IF iter >  n_iter
    #
    iter_used = min(iter-1,n_iter)
    #
    # ipos is the position in which results from the present iteraction
    # are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    #
    ipos = iter - 1 - floor(Int64, (iter-2)/n_iter)*n_iter

    println("mix_broyden: ipos      = ", ipos)
    println("mix_broyden: iter_used = ", iter_used)

    @views deltaout[:] = deltaout[:] - deltain[:]

    # Preconditioning
    CellVolume = pw.CellVolume
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    # This should be done over smooth components only

    rs = ( 3.0*CellVolume / (4*pi) / Nelectrons )^(1/3)
    agg0 = (12/pi)^(2/3) / rs
    for ig in 1:Ng
        deltaout[ig] = deltaout[ig] * G2[ig] / ( G2[ig] + agg0 )
    end


    if iter > 1
        @views df[:,ipos] = deltaout[:] - df[:,ipos]
        @views dv[:,ipos] = deltain[:]  - dv[:,ipos]
        @views nrm = sqrt( norm(df[:,ipos])^2 )
        @views df[:,ipos] = df[:,ipos]/nrm
        @views dv[:,ipos] = dv[:,ipos]/nrm
    end

    beta = zeros(maxter,maxter)

    for i in 1:iter_used
        for j in i+1:iter_used
            beta[i,j] = wg[i] * wg[j] * real(dot(df[:,j],df[:,i]))
            beta[j,i] = beta[i,j]
        end
        beta[i,i] = wg0^2 + wg[i]^2
    end

    beta_inv = inv(beta[1:iter_used,1:iter_used])

    @views beta[1:iter_used,1:iter_used] = beta_inv[:,:]

    work = zeros(iter_used)
    for i in 1:iter_used
        work[i] = real(dot(df[:,i], deltaout))
    end
    
    @views deltain[:] = deltain[:] + alphamix*deltaout[:]

    for i in 1:iter_used
        gammamix = 0.0
        for j in 1:iter_used
            gammamix = gammamix + beta[j,i] * wg[j] * work[j]
        end
        @views deltain[:] = deltain[:] - wg[i]*gammamix*( alphamix*df[:,i] + dv[:,i] )
    end

    inext = iter - floor(Int64, (iter - 1)/n_iter)*n_iter

    println("mix_broyden: inext   = ", inext)

    @views df[:,inext] = deltaout[:]
    @views dv[:,inext] = deltainsave[:]

    return

end