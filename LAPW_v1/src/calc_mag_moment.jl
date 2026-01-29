function calc_mag_moment!( atoms, pw, mt_vars, cfunir, elec_chgst, magmt, magir )

    ncmag = false # XXX Hardcoded

    mommt = elec_chgst.mommt
    mommttot = elec_chgst.mommttot
    momir = elec_chgst.momir
    momtot = elec_chgst.momtot

    wrmt = mt_vars.wrmt
    CellVolume = pw.CellVolume
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nrmt = mt_vars.nrmt
    nrmtmax = maximum(nrmt)

    ndmag = size(magir, 2) # take from here
    Npoints = size(magir, 1)

    fr = zeros(Float64, nrmtmax)
    y00 = 0.28209479177387814347

    # find the muffin-tin moments
    fill!(mommttot, 0.0)
    for idm in 1:ndmag
        for ia in 1:Natoms
            isp = atm2species[ia]
            nr = nrmt[isp]
            # extract the l=m=0 component from the muffin-tin magnetization
            @views rf_mt_lm!(1, isp, mt_vars, magmt[ia][:,idm], fr[1:nr])
            # integrate to the muffin-tin radius
            t1 = dot( wrmt[isp], fr[1:nr])
            mommt[idm,idm] = 4π * y00 * t1
            mommttot[idm] += mommt[idm,ia]
        end
    end
    # find the interstitial moments
    for idm in 1:ndmag
        t1 = dot( magir[:,idm], cfunir )
        momir[idm] = t1*CellVolume/Npoints
    end
    @. momtot[:] = mommttot[:] + momir[:]
    # total moment magnitude
    if ncmag
        elec_chgst.momtotm = sqrt(momtot[1]^2 + momtot[2]^2 + momtot[3]^2)
    else 
        elec_chgst.momtotm = abs(momtot[1])
    end
    return
end