function rhomag!(
    atoms, atsp_vars, pw, sym_vars, mt_vars, apwlo_vars, core_states, elec_chgst,
    cfunir,
    rhomt, rhoir;
    magmt = nothing, magir = nothing
)

    Natoms = atoms.Natoms
    spinpol = elec_chgst.spinpol

    if isnothing(magmt)
        ndmag = 0
    else
        ndmag = size(magir, 2)
    end

    # rho, muffin tin
    for ia in 1:Natoms
        fill!(rhomt[ia], 0.0)
    end
    # interstitial
    fill!(rhoir, 0.0)
    if spinpol
        # magnetization, muffin tin
        for ia in 1:Natoms
            fill!(magmt[ia], 0.0)
        end
        # interstitial
        fill!(magir, 0.0)
    end

    for ik in 1:pw.gvecw.kpoints.Nkpt

        evecfv = deserialize("evecs_1st_ik_$(ik).dat");
        if elec_chgst.tevecsv
            evecsv = deserialize("evecs_2nd_ik_$(ik).dat");
        else
            evecsv = nothing
        end
        apwalm = calc_match_coeffs(ik, atoms, pw, mt_vars, apwlo_vars);

        rhomagk!(
            ik, atoms, pw, mt_vars, apwlo_vars, elec_chgst,
            apwalm, evecfv, evecsv,
            rhomt, rhoir;
            magmt=magmt, magir=magir
        )

    end

    rhomagsh!(atoms, mt_vars, rhomt; magmt=magmt)

    # Symmetrize, using coarse grid
    symrfmt!(atoms, mt_vars, sym_vars, rhomt; coarse=true)
    symrfir!(pw, sym_vars, rhoir)

    rf_mt_c_to_f!(atoms, atsp_vars, mt_vars, rhomt)

    symrvfmt!(atoms, sym_vars, mt_vars, magmt; tspin = true, coarse = true)
    symrvfir!(pw, sym_vars, magir; tspin = true, tnc = false)

    for i in 1:ndmag
        @views rf_mt_c_to_f!(atoms, atsp_vars, mt_vars, magmt[:,i])
    end

    rhocore!(atoms, mt_vars, elec_chgst, core_states, rhomt; magmt=magmt)

    calc_charge!(atoms, pw, mt_vars, elec_chgst, rhomt, rhoir, cfunir)

    rhonorm!(atoms, pw, mt_vars, elec_chgst, rhomt, rhoir)

    if spinpol
        calc_mag_moment!( atoms, pw, mt_vars, cfunir, elec_chgst, magmt, magir )
    end

    return
end