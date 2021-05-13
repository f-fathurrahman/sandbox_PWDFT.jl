# Generates the coarse and fine radial meshes for each atomic species in the
# crystal. Also determines which points are in the inner part of the
# muffin-tin using the value of {\tt fracinr}.
function genrmesh!(atm_vars, atsp_vars, mt_vars)

    Nspecies = atm_vars.Nspecies

    nrsp = atsp_vars.nrsp
    rminsp = atsp_vars.rminsp
    rmaxsp = atsp_vars.rmaxsp
    
    rmt = mt_vars.rmt
    nrmt = mt_vars.nrmt

    # estimate the number of radial mesh points to infinity
    for is in 1:Nspecies
        # logarithmic mesh
        t1 = log(rmaxsp[is]/rminsp[is])/log(rmt[is]/rminsp[is])
        t2 = (nrmt[is] - 1)*t1
        nrsp[is] = round(Int64,t2) + 1 # XXX compare with nint
    end

    lmaxo = mt_vars.lmaxo

    # The following actually allocate memory
    for isp in 1:Nspecies
        atsp_vars.rsp[isp] = zeros(Float64,nrsp[isp])
        atsp_vars.rhosp[isp] = zeros(Float64, nrsp[isp])
        atsp_vars.vrsp[isp] = zeros(Float64, nrsp[isp])
    end
    
    for isp in 1:Nspecies
        mt_vars.rlmt[isp] = OffsetArray(
            zeros(Float64, nrmt[isp], 2*(lmaxo+2)),
            1:nrmt[isp], -lmaxo-1:lmaxo+2
        )
        mt_vars.wrmt[isp] = zeros(Float64, nrmt[isp])
        mt_vars.wprmt[isp] = zeros(Float64, 4, nrmt[isp])
    end

    rsp = atsp_vars.rsp
    rlmt = mt_vars.rlmt
    wrmt = mt_vars.wrmt
    wprmt = mt_vars.wprmt

    for isp in 1:Nspecies
        t1 = 1.0/(nrmt[isp] - 1)
        # logarithmic mesh
        t2 = log(rmt[isp]/rminsp[isp])
        for ir in 1:nrsp[isp]
            rsp[isp][ir] = rminsp[isp]*exp( (ir-1) * t1 * t2)
        end
        # calculate r^l on the fine radial mesh
        nr = nrmt[isp]
        for ir in 1:nr
            rlmt[isp][ir,-1] = 1.0/rsp[isp][ir]
            rlmt[isp][ir,0] = 1.0
            rlmt[isp][ir,1] = rsp[isp][ir]
        end
        #
        for l in range(-2,stop=-lmaxo-1,step=-1)
            for ir in 1:nr
                rlmt[isp][ir,l] = rlmt[isp][ir,l+1]/rsp[isp][ir]
            end
        end
        #
        for l in 2:lmaxo+2            
            for ir in 1:nr
                rlmt[isp][ir,l] = rlmt[isp][ir,l-1] * rsp[isp][ir]
            end
        end
        # determine the weights for spline integration on the fine radial mesh
        wsplint!(nr, rsp[isp], wrmt[isp])
        # multiply by r^2
        for ir in 1:nr
            wrmt[isp][ir] = wrmt[isp][ir]*rlmt[isp][ir,2]
        end
        # determine the weights for partial integration on fine radial mesh
        wsplintp!(nr, rsp[isp], wprmt[isp])
    end


    # determine the fraction of the muffin-tin radius which defines the inner part
    if mt_vars.fracinr < 0.0
        # be explicit about conversion to Float64 (not really needed actually for Julia)
        mt_vars.fracinr = sqrt( Float64(lmmaxi) / Float64(lmmaxo) )
    end

    fracinr = mt_vars.fracinr
    nrcmt = mt_vars.nrcmt
    nrmti = mt_vars.nrmti
    nrcmti = mt_vars.nrcmti
    lradstp = mt_vars.lradstp
    println("nrcmt    = ", nrcmt[1:Nspecies])

    # set up the coarse radial meshes and find the inner part of the muffin-tin
    # where rho is calculated with lmaxi
    for isp in 1:Nspecies
        mt_vars.rcmt[isp] = zeros(Float64, nrcmt[isp])
        mt_vars.rlcmt[isp] = OffsetArray(
            zeros(Float64, nrcmt[isp], 2*(lmaxo+2)),
            1:nrcmt[isp], -lmaxo-1:lmaxo+2
        )
        mt_vars.wrcmt[isp] = zeros(Float64, nrcmt[isp])
        mt_vars.wprcmt[isp] = zeros(Float64, 4, nrcmt[isp])
    end

    rcmt = mt_vars.rcmt
    rlcmt = mt_vars.rlcmt
    wrcmt = mt_vars.wrcmt
    wprcmt = mt_vars.wprcmt

    for isp in 1:Nspecies
        t1 = fracinr*rmt[isp]
        nrmti[isp] = 1
        nrcmti[isp] = 1
        irc = 0
        for ir in range(1,stop=nrmt[isp],step=lradstp)
            irc = irc + 1
            rcmt[isp][irc] = rsp[isp][ir]
            if rsp[isp][ir] < t1
                nrmti[isp] = ir
                nrcmti[isp] = irc
            end
        end
        # store r^l on the coarse radial mesh
        for l in range(-lmaxo-1, stop=lmaxo+2)
            irc = 0
            for ir in range(1,stop=nrmt[isp],step=lradstp)
                irc = irc + 1
                rlcmt[isp][irc,l] = rlmt[isp][ir,l]
            end
        end
        # determine the weights for spline integration on the coarse radial mesh
        nrc = nrcmt[isp]
        wsplint!(nrc, rcmt[isp], wrcmt[isp])
        # multiply by r^2
        #il1 = lmaxo2idx(2,lmaxo)
        for ir in 1:nrc
            wrcmt[isp][ir] = wrcmt[isp][ir] * rlcmt[isp][ir,2]
        end
        # determine the weights for partial integration on coarse radial mesh
        wsplintp!(nrc, rcmt[isp], wprcmt[isp])
    end

    return
end 


