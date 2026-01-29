# originally rfinp
function rf_inner_prod(rfmt1, rfir1, rfmt2, rfir2)
    # interstitial contribution
    res = 0.0
    for ip in 1:Npoints
        res += rfir1[ip]*rfir2[ip]*cfunir[ip]
    end
    res = res*CellVolume/Npoints

    # muffin-tin contribution
    for ia in 1:Natoms
        isp = atm2species[ia]
        res += rf_mt_inner_prod( nrmt[isp], nrmti[isp], wrmt[isp], rfmt1[ia], rfmt2[ia] )
    end
    return res
end
