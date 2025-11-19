function calc_charge!(atoms, pw, mt_vars, elec_chgst, rhomt, rhoir, cfunir)
    
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume

    # find the muffin-tin charges
    calc_chargemt!(atoms, mt_vars, elec_chgst, rhomt)
    # find the interstitial charge
    t1 = dot(rhoir, cfunir)
    elec_chgst.chgir = t1*CellVolume/Npoints
    # total calculated charge
    elec_chgst.chgcalc = elec_chgst.chgmttot + elec_chgst.chgir
    return
end

function calc_chargemt!(atoms, mt_vars, elec_chgst, rhomt)
    
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    wrmt = mt_vars.wrmt
    nrmtmax = maximum(nrmt)

    chgmt = elec_chgst.chgmt

    # automatic arrays
    fr = zeros(Float64, nrmtmax)

    y00 = 0.28209479177387814347

    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        # extract the l=m=0 component from the muffin-tin density
        @views rf_mt_lm!(1, isp, mt_vars, rhomt[ia], fr[1:nr])
        # integrate to the muffin-tin radius
        t1 = dot( wrmt[isp], fr[1:nr] )
        chgmt[ia] = 4π * y00 * t1
    end
    elec_chgst.chgmttot = sum(chgmt)
    return
end

