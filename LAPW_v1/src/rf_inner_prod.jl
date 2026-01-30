# originally rfinp
function rf_inner_prod(atoms, pw, mt_vars, cfunir, rfmt1, rfir1, rfmt2, rfir2)

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume

    # interstitial contribution
    res = 0.0
    for ip in 1:Npoints
        res += rfir1[ip]*rfir2[ip]*cfunir[ip]
    end
    res = res*CellVolume/Npoints

    # muffin-tin contribution
    for ia in 1:Natoms
        isp = atm2species[ia]
        res += rf_mt_inner_prod( isp, mt_vars, rfmt1[ia], rfmt2[ia] )
    end
    return res
end
