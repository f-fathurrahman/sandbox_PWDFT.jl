function loop_calc_qradG(
    pw::PWGrid,
    pspots::Vector{PsPot_UPF}
)
    CellVolume = pw.CellVolume
    Nspecies = length(pspots)
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

    #ndm = max( upf(:)%kkbeta )
    #nqxq = INT( ( (SQRT(ecutrho) + qnorm) / dq + 4) * cell_factor )
    # FIXME: Using floor instead of int?
    nqxq = floor(Int64, sqrt(2*ecutrho)/dq + 4) # convert to Ry
    println("nqxq = ", nqxq)


    for isp in 1:Nspecies

        psp = pspots[isp] # shorthand

        # skip if not using ultrasoft
        if !pspots[isp].is_ultrasoft
            continue
        end

        println("\n-----------------------")
        @printf("Species = %3d nqlc = %3d nbeta = %3d\n", isp, psp.nqlc, psp.Nproj)
        println("kkbeta  = ", psp.kkbeta)
        println("size qfuncl = ", size(psp.qfuncl))
        println("lmax = ", psp.lmax)

        for l in 0:psp.nqlc-1
            println("l = ", l)
            # note that l is the true (combined) angular momentum
            # and that the arrays have dimensions 0..l (no more 1..l+1)
            #for iq in 1:nqxq
            for iq in 1:1
                q = (iq - 1) * dq
                # here we compute the spherical bessel function for each q_i
                #CALL sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, q, l, besr)
                #
                for nb in 1:psp.Nproj
                    # the Q are symmetric with respect to indices
                    for mb in nb:psp.Nproj
                        ijv = round(Int64, mb*(mb - 1)/2) + nb
                        lnb = psp.proj_l[nb]
                        lmb = psp.proj_l[mb]
                        cond1 = l >= abs( lnb - lmb )
                        cond2 = l <=      lnb + lmb
                        cond3 = mod(l + lnb + lmb, 2) == 0
                        @printf("nb,lnb=(%3d,%3d) mb,lnb=(%3d,%3d) ijv=%3d", nb, lnb, mb, lmb, ijv)
                        if cond1 & cond2 & cond3
                            println(" qrad is evaluated")
                            #for ir = 1, upf(nt)%kkbeta
                            #    aux[ir] = besr(ir) * upf(nt)%qfuncl(ir,ijv,l)
                            #end
                            #
                            # and then we integrate with all the Q functions
                            #CALL simpson( upf(nt)%kkbeta, aux, rgrid(nt)%rab, qrad(iq,ijv,l+1,nt) )
                        else
                            println()
                        end
                    end
                end # igl
            end # l
        end
        # qrad(:,:,:,nt) = qrad(:,:,:,nt)*prefr
        # CALL mp_sum( qrad(:,:,:,nt), intra_bgrp_comm )
    end
    return
end
