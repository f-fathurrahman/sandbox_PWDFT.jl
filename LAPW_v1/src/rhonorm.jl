function rhonorm!(atoms, pw, mt_vars, elec_chgst, rhomt, rhoir)
    
    epschg = 0.001
    y00 = 0.28209479177387814347

    chgcalc = elec_chgst.chgcalc
    chgtot = elec_chgst.chgtot
    chgmt = elec_chgst.chgmt
    CellVolume = pw.CellVolume # XXX: add CellVolume to Atoms?
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    rmt = mt_vars.rmt

    # check error in total charge
    t1 = chgcalc/chgtot - 1.0
    println("diff charge = $t1")
    if abs(t1) > epschg
        println("WARNING: Total charge density is incorrect: calculated: $chgcalc, required: $chgtot")
    end

    # error in average density
    t1 = (chgtot - chgcalc)/CellVolume
    
    # add the constant difference to the density
    t2 = t1/y00
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        i = 1
        for ir in 1:nri
            rhomt[ia][i] += t2
            i += lmmaxi
        end
        for ir in (nri+1):nr
            rhomt[ia][i] += t2
            i += lmmaxo
        end
    end
    rhoir .+= t1
    # add the difference to the charges
    for isp in 1:Nspecies
        t2 = t1*(4π/3)*rmt[isp]^3
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            chgmt[ia] += t2
            elec_chgst.chgmttot += t2
        end
    end
    elec_chgst.chgir = elec_chgst.chgtot - elec_chgst.chgmttot
    return
end
