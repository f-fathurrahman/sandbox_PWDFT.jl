function calc_match_coeffs!(ik, atoms, pw, mt_vars, apwlo_vars, apwalm)
    # XXX: Need to test in case of omax > 1
    # apword must be at least 2, however it is quite difficult to figure out
    # apwe0 for this.

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Nspecies = atoms.Nspecies

    apwordmax = apwlo_vars.apwordmax # accross all species

    lmaxapw = mt_vars.lmaxapw
    lmmaxapw = mt_vars.lmmaxapw
    nrmt = mt_vars.nrmt
    rmt = mt_vars.rmt
    idxlm = mt_vars.idxlm

    apwfr = apwlo_vars.apwfr

    Ngk = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g
    kpoints = pw.gvecw.kpoints
    vgpc = zeros(Float64, 3, Ngk)
    gpc = zeros(Float64, Ngk)
    @views kvec = kpoints.k[:,ik]
    for igk in 1:Ngk
        ig = idx_gw2g[ik][igk]
        vgpc[1:3,igk] .= pw.gvec.G[1:3,ig] .+ kvec[1:3]
        gpc[igk] = norm(vgpc[1:3,igk])
    end
    sfacgp = zeros(ComplexF64, Ngk, Natoms)
    gensfacgp!(atoms, Ngk, vgpc, sfacgp)

    a = zeros(ComplexF64, apwordmax, apwordmax)

    djl = OffsetArray(
        zeros(Float64, lmaxapw+1, apwordmax, Ngk),
        0:lmaxapw, 1:apwordmax, 1:Ngk
    )
    ylmgp = zeros(ComplexF64, lmmaxapw, Ngk)

    if apwordmax > 1
        b = zeros(ComplexF64, apwordmax, Ngk*(2*lmaxapw+1))
    end
  
    # compute the spherical harmonics of the G+p-vectors
    for igp in 1:Ngk
        @views genylmv!( lmaxapw, vgpc[:,igp], ylmgp[:,igp] )
    end

    t0 = 4π/sqrt(pw.CellVolume)
  
    # loop over species
    for isp in 1:Nspecies
        nr = nrmt[isp]
        # maximum APW order for this species
        omax = maximum( apwlo_vars.apword[isp] )
        # special case of omax=1
        if omax==1 
            for igp in 1:Ngk
                t1 = gpc[igp]*rmt[isp]
                @views sbessel!(lmaxapw, t1, djl[:,1,igp] )
            end
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                for l in 0:lmaxapw
                    z1 = ( t0/apwfr[ia][l][1][nr,1] ) * im^(l) #zil(l)
                    for igp in 1:Ngk
                        z2 = djl[l,1,igp] * z1 * sfacgp[igp,ia]
                        for m in -l:l
                            lm = idxlm[l,m]
                            apwalm[ia][igp,1,lm] = z2*conj(ylmgp[lm,igp])
                        end 
                    end 
                end 
            end 
            continue # CYCLE  # next species
        end # if omax == 1 
        # starting point on radial mesh for fitting polynomial of order npapw
        ir = nr-npapw+1
        # evaluate the spherical Bessel function derivatives for all G+p-vectors
        for igp in 1:Ngk
            t1 = gpc[igp]*rmt[isp]
            for io in 1:omax
                @views sbesseldm!( io-1, lmaxapw, t1, djl[:,io,igp] )
            end 
            t1 = 1.0
            for io in 2:omax
                t1 = t1*gpc[igp]
                djl[:,io,igp] = t1*djl[:,io,igp]
            end 
        end 
        # loop over atoms
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            # begin loop over l
            for l in 0:lmaxapw
                z1 = t0*im^l  #zil[l]
                # set up matrix of derivatives
                for jo in 1:apword[isp][l], io in 1:apword[isp][l]
                    a[io,jo] = polynm(io-1, npapw, rsp[isp][ir], apwfr[ia][l][jo][ir,1], rmt[isp])
                end 
                # set up target vectors
                i = 0
                for igp in 1:Ngk
                    z2 = z1*sfacgp[igp,ia]
                    for m in -l:l
                        lm = idxlm(l,m)
                        i = i + 1
                        z3 = z2*conjg(ylmgp[lm,igp])
                        for io in 1:apword[isp][l]
                            b[io,i] = djl[l,io,igp]*z3
                        end 
                    end 
                end 
                # solve the general complex linear systems
                #CALL zgesv(apword(l,is),i,a,apwordmax,ipiv,b,apwordmax,info)
                idx1 = 1:apword[isp][l]
                b[idx1,idx1] = a[idx1,idx1] \ b[idx1,idx1]
                #
                i = 0
                for igp in 1:Ngk
                    for m in -l:l
                        lm = idxlm[l,m]
                        i = i + 1
                        for io in 1:apword[isp][l]
                            apwalm[ia][igp,io,lm] = b[io,i]
                        end
                    end 
                end 
            end  # end loop over l 
        end 
    end # end loops over atoms and species

    return
end

