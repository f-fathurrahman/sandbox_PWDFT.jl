function calc_qradG(atoms, pw, pspots)
    CellVolume = pw.CellVolume
    Nspecies = atoms.Nspecies
    prefr = 4π/CellVolume
    ecutrho = pw.ecutrho

    ndm = 0
    for isp in 1:Nspecies
        println("kkbeta = ", pspots[isp].kkbeta)
        if ndm < pspots[isp].kkbeta
            ndm = pspots[isp].kkbeta
        end
    end
    println("ndm = ", ndm)

    qnorm = 0.0 # XXX HARDCODED, no k-points norm of (q + k) ?
    dq = 0.01 # XXX HARDCODED
    cell_factor = 1.0 # hardcoded
    nqxq = floor(Int64, sqrt(2*ecutrho)/dq + 4) # factor of 2 in 2*ecutrho (convert to Ry)
    println("nqxq = ", nqxq)


    besr = zeros(Float64, ndm)
    aux = zeros(Float64, ndm)
    # Radial Fourier transform of
    qradG = Array{Array{Float64,3},1}(undef,Nspecies)
    for isp in 1:Nspecies
        Nproj = pspots[isp].Nproj
        Nn2 = round(Int64, Nproj*(Nproj+1)/2)
        # Determine lmaxq
        lmaxkb = -1
        for i in 1:pspots[isp].Nproj
            lmaxkb = max(lmaxkb, pspots[isp].proj_l[i])
        end
        lmaxq = 2*lmaxkb + 1
        qradG[isp] = zeros(Float64, nqxq, Nn2, lmaxq)
    end
    # third dimension can be changed to (lmaxkb[isp] + 1)
    # at least lmaxq=1


    for isp in 1:Nspecies
        #
        psp = pspots[isp] # shorthand
        #
        # skip if not using ultrasoft
        if !pspots[isp].is_ultrasoft
            continue
        end

        for l in 0:psp.nqlc-1
            # note that l is the true (combined) angular momentum
            # and that the arrays have dimensions 0..l (no more 1..l+1)
            for iq in 1:nqxq
                q = (iq - 1) * dq
                # here we compute the spherical bessel function for each q_i
                for ir in 1:psp.kkbeta
                    besr[ir] = sphericalbesselj(l, q*psp.r[ir])
                end
                #
                for nb in 1:psp.Nproj
                    # the Q are symmetric w.r.t indices
                    for mb in nb:psp.Nproj
                        ijv = round(Int64, mb*(mb - 1)/2) + nb
                        lnb = psp.proj_l[nb]
                        lmb = psp.proj_l[mb]
                        cond1 = l >= abs( lnb - lmb )
                        cond2 = l <=      lnb + lmb
                        cond3 = mod(l + lnb + lmb, 2) == 0
                        # XXX: Check qfuncl for this conditions?
                        if cond1 & cond2 & cond3
                            for ir in 1:psp.kkbeta
                                aux[ir] = besr[ir] * psp.qfuncl[ir,ijv,l+1]
                            end
                            # and then we integrate with all the Q functions
                            qradG[isp][iq,ijv,l+1] = PWDFT.integ_simpson( psp.kkbeta, aux, psp.rab )
                        end
                    end
                end # igl
            end # l
        end
        qradG[isp][:,:,:] = qradG[isp][:,:,:]*prefr
    end
    return qradG
end



