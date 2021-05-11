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

    atsp_vars.rsp = zeros(Float64, nrspmax, Nspecies)
    mt_vars.rlmt = zeros(Float64, nrmtmax, 2*(lmaxo+2), Nspecies)
    mt_vars.wrmt = zeros(Float64, nrmtmax, Nspecies)
    mt_vars.wprmt = zeros(Float64, 4, nrmtmax, Nspecies)

    rsp = atsp_vars.rsp
    rlmt = mt_vars.rlmt
    wrmt = mt_vars.wrmt
    wprmt = mt_vars.wprmt

    #println("2*(lmaxo+2) = ", 2*(lmaxo+2))
    #println("size(rlmt) = ", size(rlmt))

    for is in 1:Nspecies
        t1 = 1.0/(nrmt[is] - 1)
        # logarithmic mesh
        t2 = log(rmt[is]/rminsp[is])
        for ir in 1:nrsp[is]
            rsp[ir,is] = rminsp[is]*exp( (ir-1) * t1 * t2)
        end
        # calculate r^l on the fine radial mesh
        nr = nrmt[is]
        # XXX We are using lmaxo2idx
        il1 = lmaxo2idx(-1,lmaxo)
        il2 = lmaxo2idx( 0,lmaxo)
        il3 = lmaxo2idx( 1,lmaxo)
        for ir in 1:nr
            rlmt[ir,il1,is] = 1.0/rsp[ir,is]
            rlmt[ir,il2,is] = 1.0
            rlmt[ir,il3,is] = rsp[ir,is]
        end
        #
        for l in range(-2,stop=-lmaxo-1,step=-1)
            il1 = lmaxo2idx(l,lmaxo)
            il2 = lmaxo2idx(l+1,lmaxo)
            for ir in 1:nr
                rlmt[ir,il1,is] = rlmt[ir,il2,is]/rsp[ir,is]
            end
        end
        #
        for l in 2:lmaxo+2            
            il1 = lmaxo2idx(l,lmaxo)
            il2 = lmaxo2idx(l-1,lmaxo)
            for ir in 1:nr
                rlmt[ir,il1,is] = rlmt[ir,il2,is] * rsp[ir,is]
            end
        end
        # determine the weights for spline integration on the fine radial mesh
        @views wsplint!(nr, rsp[:,is], wrmt[:,is])
        # multiply by r^2
        il2 = lmaxo2idx(2,lmaxo)
        for ir in 1:nr
            wrmt[ir,is] = wrmt[ir,is]*rlmt[ir,il2,is]
        end
        # determine the weights for partial integration on fine radial mesh
        @views wsplintp!(nr, rsp[:,is], wprmt[:,:,is])
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
    mt_vars.rlcmt = zeros(Float64, nrcmtmax, 2*(lmaxo+2), Nspecies)
    mt_vars.wrcmt = zeros(Float64, nrcmtmax, Nspecies)
    mt_vars.wprcmt = zeros(Float64, 4, nrcmtmax, Nspecies)

    rcmt = mt_vars.rcmt
    rlcmt = mt_vars.rlcmt
    wrcmt = mt_vars.wrcmt
    wprcmt = mt_vars.wprcmt

    for is in 1:Nspecies
        t1 = fracinr*rmt[is]
        nrmti[is] = 1
        nrcmti[is] = 1
        irc = 0
        for ir in range(1,stop=nrmt[is],step=lradstp)
            irc = irc + 1
            rcmt[irc,is] = rsp[ir,is]
            if rsp[ir,is] < t1
                nrmti[is] = ir
                nrcmti[is] = irc
            end
        end
        # store r^l on the coarse radial mesh
        for l in range(-lmaxo-1, stop=lmaxo+2)
            irc = 0
            il1 = lmaxo2idx(l,lmaxo)
            for ir in range(1,stop=nrmt[is],step=lradstp)
                irc = irc + 1
                rlcmt[irc,il1,is] = rlmt[ir,il1,is]
            end
        end
        # determine the weights for spline integration on the coarse radial mesh
        nrc = nrcmt[is]
        @views wsplint!(nrc, rcmt[:,is], wrcmt[:,is])
        # multiply by r^2
        il1 = lmaxo2idx(2,lmaxo)
        for ir in 1:nrc
            wrcmt[ir,is] = wrcmt[ir,is] * rlcmt[ir,il1,is]
        end
        # determine the weights for partial integration on coarse radial mesh
        @views wsplintp!(nrc, rcmt[:,is], wprcmt[:,:,is])
    end

    return
end 


