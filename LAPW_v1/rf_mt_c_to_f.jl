function rf_mt_c_to_f!(
    atoms, mt_vars, rfmt
)
    # USE m_atoms, ONLY: natmtot, rsp, idxis
    # USE m_muffin_tins, ONLY: npmtmax, nrcmtmax, nrmtmax, npmt, rcmt, npcmti, nrcmt, nrcmti, &
    #                  npmti, lmmaxo, lmmaxi, nrmti, nrmt, lradstp, rlmt
    # IMPLICIT NONE 
    # ! arguments
    # REAL(8), intent(inout) :: rfmt(npmtmax,natmtot)
    # ! local variables
    # INTEGER :: is,ias,lm
    # INTEGER :: nr,nri,nro
    # INTEGER :: iro,ir,npi,i
    # INTEGER :: nrc,nrci,nrco
    # INTEGER :: irco,irc,npci
    # ! ALLOCATABLE arrays
    # REAL(8), ALLOCATABLE :: fi(:),fo(:),rfmt1(:)

    if mt_vars.lradstp == 1
        return
    end

    idx2species = atoms.idx2species
    Natoms = atoms.Natoms

    nrcmt = mt_vars.nrcmt
    nrmt = mt_vars.nrmt
    nrcmti = mt_vars.nrcmti
    npmt = mt_vars.npmt
    npcmti = mt_vars.npcmti
    lmmaxo = mt_vars.lmmaxo
    lmmaxi = mt_vars.lmmaxi

    nrcmtmax = maximum(nrcmt)
    nrmtmax = maximum(nrmt)
    npmtmax = maximum(npmtmax) 

    fi = zeros(Float64, nrcmtmax)
    fo = zeros(Float64, nrmtmax)
    rfmt1 = zeros(Float64, npmtmax)
    
    for ia in 1:Natoms
        isp = idx2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        nro = nr - nri
        iro = nri + 1
        npi = npmti[isp]
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        nrco = nrc - nrci
        irco = nrci + 1
        npci = npcmti[is]
        # interpolate up to lmaxi over entire muffin-tin
        for lm in 1:lmmaxi
            i = lm
            for irc in 1:nrci
                fi[irc] = rfmt[isp][i]
                i = i + lmmaxi
            end
            for irc in irco:nrc
              fi[irc] = rfmt[isp][i]
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
                fi[irc] = rfmt[isp][i]
                i = i + lmmaxo
            end
            idxc = irco:nrco-1
            idx = iro:nro-1
            @views rf_interp(
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

    end 
    #DEALLOCATE(fi,fo,rfmt1)

    return
end
