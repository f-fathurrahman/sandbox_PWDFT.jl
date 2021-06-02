function rf_mt_c_to_f!(
    atoms, atsp_vars, mt_vars, rfmt
)
    
    # No need to interpolate if not using coarse/fine radial grid
    if mt_vars.lradstp == 1
        return
    end

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms

    nrcmt = mt_vars.nrcmt
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    nrcmti = mt_vars.nrcmti
    npmt = mt_vars.npmt
    npmti = mt_vars.npmti
    npcmti = mt_vars.npcmti
    lmmaxo = mt_vars.lmmaxo
    lmmaxi = mt_vars.lmmaxi
    rcmt = mt_vars.rcmt
    rlmt = mt_vars.rlmt

    rsp = atsp_vars.rsp

    nrcmtmax = maximum(nrcmt)
    nrmtmax = maximum(nrmt)
    npmtmax = maximum(npmt) 

    fi = zeros(Float64, nrcmtmax)
    fo = zeros(Float64, nrmtmax)
    rfmt1 = zeros(Float64, npmtmax)
    
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        nro = nr - nri
        iro = nri + 1
        npi = npmti[isp]
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        nrco = nrc - nrci
        irco = nrci + 1
        npci = npcmti[isp]
        # interpolate up to lmaxi over entire muffin-tin
        for lm in 1:lmmaxi
            i = lm
            for irc in 1:nrci
                fi[irc] = rfmt[ia][i]
                i = i + lmmaxi
            end
            for irc in irco:nrc
              fi[irc] = rfmt[ia][i]
              i = i + lmmaxo
            end
            @views rf_interp!(
                nrc, rcmt[isp][1:nrc], fi[1:nrc],
                nr, rlmt[isp][1:nr,1], fo[1:nr]
            )
            i = lm
            for ir in 1:nri
                rfmt1[i] = fo[ir]
                i = i + lmmaxi
            end
            for ir in iro:nr
              rfmt1[i] = fo[ir]
              i = i + lmmaxo
            end
        end
        # interpolate up to lmaxo on outer part of muffin-tin
        for lm in lmmaxi+1:lmmaxo
            i = npci + lm
            for irc in irco:nrc
                fi[irc] = rfmt[ia][i]
                i = i + lmmaxo
            end
            idxc = irco:irco+nrco-1
            idx = iro:iro+nro-1
            @views rf_interp!(
                nrco, rcmt[isp][idxc], fi[idxc],
                nro, rsp[isp][idx], fo[idx]
            )
            i = npi + lm
            for ir in iro:nr
              rfmt1[i] = fo[ir]
              i = i + lmmaxo
            end
        end 
        # CALL dcopy(npmt(is), rfmt1,1, rfmt(:,ias), 1)
        for ip in 1:npmt[isp]
            rfmt[ia][ip] = rfmt1[ip]
        end
    end 
    return
end
