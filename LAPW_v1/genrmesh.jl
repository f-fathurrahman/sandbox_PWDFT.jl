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
    atsp_vars.nrspmax = 1
    for is in 1:Nspecies
        #println("rminsp = ", rminsp[is])
        #println("rmaxsp = ", rmaxsp[is])
        #println("nrmt   = ", nrmt[is])
        # logarithmic mesh
        t1 = log(rmaxsp[is]/rminsp[is])/log(rmt[is]/rminsp[is])
        t2 = (nrmt[is] - 1)*t1
        #println("t1 = ", t1)
        #println("t2 = ", t2)
        nrsp[is] = round(Int64,t2) + 1 # XXX compare with nint
        atsp_vars.nrspmax = max(atsp_vars.nrspmax, nrsp[is])
    end

    nrspmax = atsp_vars.nrspmax
    lmaxo = mt_vars.lmaxo
    nrmtmax = mt_vars.nrmtmax
    #println("lmaxo = ", lmaxo)
    #println("nrmtmax = ", nrmtmax)

    # The following actually allocate memory
    atsp_vars.rsp = Vector{Vector{Float64}}(undef,Nspecies)
    for isp in 1:Nspecies
        atsp_vars.rsp[isp] = zeros(Float64,nrsp[isp])
    end
    
    mt_vars.rlmt = OffsetArray(
        zeros(Float64, nrmtmax, 2*(lmaxo+2), Nspecies),
        1:nrmtmax, -lmaxo-1:lmaxo+2, 1:Nspecies
    )

    mt_vars.wrmt = zeros(Float64, nrmtmax, Nspecies)
    mt_vars.wprmt = zeros(Float64, 4, nrmtmax, Nspecies)

    rsp = atsp_vars.rsp
    rlmt = mt_vars.rlmt
    wrmt = mt_vars.wrmt
    wprmt = mt_vars.wprmt

    #println("2*(lmaxo+2) = ", 2*(lmaxo+2))
    #println("size(rlmt) = ", size(rlmt))

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
            rlmt[ir,-1,isp] = 1.0/rsp[isp][ir]
            rlmt[ir,0,isp] = 1.0
            rlmt[ir,1,isp] = rsp[isp][ir]
        end
        #
        for l in range(-2,stop=-lmaxo-1,step=-1)
            #il1 = lmaxo2idx(l,lmaxo)
            #il2 = lmaxo2idx(l+1,lmaxo)
            for ir in 1:nr
                rlmt[ir,l,isp] = rlmt[ir,l+1,isp]/rsp[isp][ir]
            end
        end
        #
        for l in 2:lmaxo+2            
            #il1 = lmaxo2idx(l,lmaxo)
            #il2 = lmaxo2idx(l-1,lmaxo)
            for ir in 1:nr
                rlmt[ir,l,isp] = rlmt[ir,l-1,isp] * rsp[isp][ir]
            end
        end
        # determine the weights for spline integration on the fine radial mesh
        @views wsplint!(nr, rsp[isp], wrmt[:,isp])
        # multiply by r^2
        #il2 = lmaxo2idx(2,lmaxo)
        for ir in 1:nr
            wrmt[ir,isp] = wrmt[ir,isp]*rlmt[ir,2,isp]
        end
        # determine the weights for partial integration on fine radial mesh
        @views wsplintp!(nr, rsp[isp], wprmt[:,:,isp])
    end


    # determine the fraction of the muffin-tin radius which defines the inner part
    if mt_vars.fracinr < 0.0
        # be explicit about conversion to Float64 (not really needed actually for Julia)
        mt_vars.fracinr = sqrt( Float64(lmmaxi) / Float64(lmmaxo) )
    end

    fracinr = mt_vars.fracinr
    nrcmtmax = mt_vars.nrcmtmax
    nrcmt = mt_vars.nrcmt
    nrmti = mt_vars.nrmti
    nrcmti = mt_vars.nrcmti
    lradstp = mt_vars.lradstp
    println("nrcmtmax = ", nrcmtmax)
    println("nrcmt    = ", nrcmt[1:Nspecies])

    # set up the coarse radial meshes and find the inner part of the muffin-tin
    # where rho is calculated with lmaxi
    mt_vars.rcmt = zeros(Float64, nrcmtmax, Nspecies)
    mt_vars.rlcmt = OffsetArray(
        zeros(Float64, nrcmtmax, 2*(lmaxo+2), Nspecies),
        1:nrcmtmax, -lmaxo-1:lmaxo+2, 1:Nspecies
    )

    mt_vars.wrcmt = zeros(Float64, nrcmtmax, Nspecies)
    mt_vars.wprcmt = zeros(Float64, 4, nrcmtmax, Nspecies)

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
            rcmt[irc,isp] = rsp[isp][ir]
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
                rlcmt[irc,l,isp] = rlmt[ir,l,isp]
            end
        end
        # determine the weights for spline integration on the coarse radial mesh
        nrc = nrcmt[isp]
        @views wsplint!(nrc, rcmt[:,isp], wrcmt[:,isp])
        # multiply by r^2
        #il1 = lmaxo2idx(2,lmaxo)
        for ir in 1:nrc
            wrcmt[ir,isp] = wrcmt[ir,isp] * rlcmt[ir,2,isp]
        end
        # determine the weights for partial integration on coarse radial mesh
        @views wsplintp!(nrc, rcmt[:,isp], wprcmt[:,:,isp])
    end

    return
end 


