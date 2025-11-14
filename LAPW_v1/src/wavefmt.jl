
function wavefmt!(
    lrstp, ia, atoms, mt_vars, apwlo_vars, ngp, apwalm_ia, evecfv_ist, wfmt
)

    atm2species = atoms.atm2species

    lmaxo = mt_vars.lmaxo
    lmaxi = mt_vars.lmaxi
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    nrmti = mt_vars.nrmti
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    npcmt = mt_vars.npcmt
    npcmti = mt_vars.npcmti
    lradstp = mt_vars.lradstp
    idxlm = mt_vars.idxlm

    apword = apwlo_vars.apword
    apwfr = apwlo_vars.apwfr
    nlorb = apwlo_vars.nlorb
    lorbl = apwlo_vars.lorbl
    idxlo = apwlo_vars.idxlo
    lofr = apwlo_vars.lofr

    isp = atm2species[ia]
    ldi = lmmaxi # in Elk this is 2 times larger because it uses real arrays
    ldo = lmmaxo
    iro = nrmti[isp] + lrstp

    if lrstp == 1
        # not using coarse vs fine grid
        # XXX Are there any use cases for this?
        # XXX Probably for testing
        nrc = nrmt[isp]
        nrci = nrmti[isp]
        npc = npmt[isp]
        npci = npmti[isp]
    elseif lrstp == lradstp
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        npc = npcmt[isp]
        npci = npcmti[isp]
    else
        error("Invalid lrstp = $(lrstp)")
    end
    nrco = nrc - nrci
    
    #println("ldi=$ldi ldo=$ldo npc = $npc lrstp=$lrstp")
    #println("nrci=$nrci nrco=$nrco iro=$iro")
    
    #
    # zero the wavefunction
    @views wfmt[1:npc] .= 0.0
    # in the original Elk code this array have leading dimension of 2 and of type real(8)
    #
    # APW functions     
    #
    lm = 0
    for l in 0:lmaxo
        for m in -l:l
            lm = lm + 1
            i = npci + lm
            for io in 1:apword[isp][l]
                z1 = BLAS.dotu(ngp, evecfv_ist, 1, apwalm_ia[:,io,lm], 1)
                #println("lm=$lm io=$io z1 = $z1")
                if l <= lmaxi
                    # call daxpy(nrci,dble(z1),apwfr(1,1,io,l,ias),lrstp,wfmt(1,lm),ldi)
                    ix = 1
                    iy = 1
                    for ii in 1:nrci
                        #println("ii=$ii ix=$ix iy=$iy")
                        wfmt[iy] += z1 * apwfr[ia][l][io][ix,1]
                        ix += lrstp
                        iy += ldi
                    end
                end
                # call daxpy(nrco,dble(z1), apwfr(iro,1,io,l,ias),lrstp, wfmt(1,i),ldo)
                ix = iro
                iy = 1
                #println("\nsize apwfr[ia][l][io] = ", size(apwfr[ia][l][io]))
                #println("size wfmt = ", size(wfmt))
                for ii in 1:nrco
                    #println("ii=$ii ix=$ix iy=$iy")
                    wfmt[iy] += z1 * apwfr[ia][l][io][ix,1]
                    ix += lrstp
                    iy += ldo
                end
                #println(" temp sum wfmt = ", sum(wfmt))
            end
        end
    end
    #
    # local-orbital functions     
    #
    for ilo in 1:nlorb[isp]
        l = lorbl[isp][ilo]
        for m in -l:l
            lm = idxlm[l,m]
            i = npci + lm
            z1 = evecfv_ist[ngp+idxlo[ia][ilo][lm]]
            if l <= lmaxi
                # daxpy(nrci,dble(z1),lofr(1,1,ilo,ias),lrstp,wfmt(1,lm),ldi)
                ix = 1
                iy = 1
                for i in 1:nrci
                    wfmt[iy] += z1 * lofr[ia][ilo][ix,1]
                    ix += lrstp
                    iy += ldi
                end
            end
            #daxpy(nrco,dble(z1),lofr(iro,1,ilo,ias),lrstp,wfmt(1,i),ldo)
            ix = iro
            iy = 1
            for i in 1:nrco
                wfmt[iy] += z1 * lofr[ia][ilo][ix,1]
                ix += lrstp
                iy += ldo
            end
        end
    end

    return
end
