function calc_energy_terms!(
    atoms, atsp_vars, core_states,
    pw, mt_vars, elec_chgst, ndmag,
    cfunir,
    rhomt, rhoir,
    vsmt,
    vclmt, vclir,
    epsxcmt, epsxcir, vxcmt, vxcir,
    bsmt, bsir, magmt, magir, bxcmt, bxcir
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume

    nstsv = elec_chgst.nstsv
    occsv = elec_chgst.occsv
    evalsv = elec_chgst.evalsv

    spinpol = elec_chgst.spinpol

    E_vxc = rf_inner_prod(atoms, pw, mt_vars, cfunir, rhomt, rhoir, vxcmt, vxcir)

    # exchange-correlation effective field energy
    E_bxc = 0.0
    if spinpol
        for idm in 1:ndmag
            # interstitial contribution
            for ip in 1:Npoints
                E_bxc += magir[ip,idm]*bxcir[ip,idm]*cfunir[ip]
            end
            E_bxc *= CellVolume/Npoints
            # muffin-tin contribution
            for ia in 1:Natoms
                isp = atm2species[ia]
                E_bxc += rf_mt_inner_prod( isp, mt_vars, magmt[ia][:,idm], bxcmt[ia] )
            end
        end
    end

#=
    momtot = elec_chgst.momtot
    # external magnetic field energies
    E_bext = 0.0
    gfacte = 2.00231930436256
    solsc = 137.035999084 # there can be additional factor
    cb = gfacte/(4.0*solsc)
    for idm in 1:ndmag
        if ncmag
            jdm = idm
        else
            jdm = 3
        end
        # energy of physical global field
        E_bext += cb*momtot[idm]*bfieldc[jdm]
    end
=#

    E_vcl = rf_inner_prod(atoms, pw, mt_vars, cfunir, rhomt, rhoir, vclmt, vclir)
    

    y00 = 0.5/sqrt(pi)
    E_mad = 0.0
    Zatoms = -PWDFT.get_Zatoms(atoms) # use minus sign
    for ia in 1:Natoms
        isp = atm2species[ia]
        E_mad += 0.5*Zatoms[ia]*(vclmt[ia][1] - atsp_vars.vcln[isp][1])*y00
    end

    # This should be called only once in the beginning of SCF
    E_nn = calc_engynn(atoms, atsp_vars, pw, mt_vars)

    # electron-nuclear interaction energy
    E_en = 2.0*(E_mad - E_nn)
    
    # Hartree energy
    E_har = 0.5*(E_vcl - E_en)

    # Coulomb energy
    E_cl = E_nn + E_en + E_har

    E_xc = rf_inner_prod(atoms, pw, mt_vars, cfunir, rhomt, rhoir, epsxcmt, epsxcir)

    E_kin_core = calc_Ekin_core(atoms, atsp_vars, core_states, mt_vars, vsmt)

    nstsp = atsp_vars.nstsp
    occcr = core_states.occcr
    evalcr = core_states.evalcr
    spcore = atsp_vars.spcore
    # sum of eigenvalues
    # core eigenvalues
    evalsum = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ist in 1:nstsp[isp]
            if spcore[isp][ist]
                evalsum += occcr[ia][ist] * evalcr[ia][ist]
            end
        end
    end

    # valence eigenvalues
    for ik in 1:Nkpt
        for ist in 1:nstsv
            evalsum += wk[ik]*occsv[ist,ik]*evalsv[ist,ik]
        end
    end

    #
    ss_mag_field = 0.0
    if spinpol
        rfmt = Vector{Vector{Float64}}(undef, Natoms)
        # remove magnetic field contribution
        for idm in 1:ndmag
            for ia in 1:Natoms
                isp = atm2species[ia]
                rfmt[ia] = zeros(Float64, mt_vars.npmt[isp])
                forward_SHT!(mt_vars, isp, bsmt[ia][:,idm], rfmt[ia], coarse=true)
            end
            rf_mt_c_to_f!(atoms, atsp_vars, mt_vars, rfmt)
            # interstitial contribution
            for ip in 1:Npoints
                ss_mag_field += magir[ip,idm]*bsir[ip,idm]*cfunir[ip]
            end
            ss_mag_field *= CellVolume/Npoints
            # muffin-tin contribution
            for ia in 1:Natoms
                isp = atm2species[ia]
                ss_mag_field += rf_mt_inner_prod( isp, mt_vars, magmt[ia][:,idm], rfmt[ia] )
            end
        end # for
    end # if

    E_kin = evalsum - E_vcl - E_vxc - ss_mag_field
    
    # Move this to elec_chgst ?
    if elec_chgst.nspinor == 2
        occmax = 1.0
    elseif elec_chgst.nspinor == 1
        occmax = 2.0
    else
        @error("Invalid nspinor value")
    end
    #
    swidth = elec_chgst.swidth
    ha_si = 4.3597447222071e-18
    kb_si = 1.380649e-23
    kboltz = kb_si/ha_si  # Boltzmann constant in Hartree/kelvin
    # This is only for Fermi-Dirac smearing
    # entropic contribution
    entrpy = 0.0
    E_TS = 0.0
    # non-zero only for the Fermi-Dirac smearing function
    ss = 0.0
    for ik in 1:Nkpt
        for ist in 1:nstsv
            f = occsv[ist,ik]/occmax
            if (f > 0.0) && (f < 1.0)
                ss += wk[ik]*(f*log(f) + (1.0-f)*log(1.0-f))
            end
        end
    end
    # entropy
    entrpy = -occmax*kboltz*ss
    # contribution to free energy
    E_TS = -swidth*entrpy/kboltz
    
    E_tot = E_kin + 0.5*E_vcl + E_mad + E_xc + E_TS

    println("E_kin = ", E_kin)
    println("E_vxc = ", E_vxc)
    println("E_vcl = ", E_vcl)
    println("E_mad = ", E_mad)
    println("E_nn = ", E_nn)
    println("E_en = ", E_en)
    println("E_har = ", E_har)
    println("E_cl = ", E_cl)
    println("E_xc = ", E_xc)
    println("E_kin_core = ", E_kin_core)
    println("entrpy = ", entrpy)
    println("E_TS = ", E_TS)
    println("E_bxc = ", E_bxc)
    println("ss_mag_field = ", ss_mag_field)
    println("evalsum = ", evalsum)
    println("E_tot = ", E_tot)

    return E_tot
end
