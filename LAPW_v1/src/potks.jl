function potks!(
    atoms, atsp_vars, mt_vars, pw, sym_vars, 
    rhomt, rhoir,
    vclmt, vclir,
    epsxcmt, epsxcir,
    vxcmt, vxcir,
    vsmt, vsir;
    magmt = nothing, magir = nothing,
    bxcir = nothing, bxcmt = nothing,
    bsmt = nothing, bsir = nothing,
    cfunir = nothing,
    bfieldc = nothing, bfcmt = nothing,
    spinpol = false, ncmag = false
)

    Natoms = atoms.Natoms

    # Solver Hartree equation (compute electrostatic Coulomb potential)
    potcoul!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir, vclmt, vclir )

    if spinpol
        potxcmt!(atoms, mt_vars, rhomt, magmt, epsxcmt, vxcmt, bxcmt)
    else
        potxcmt!(atoms, mt_vars, rhomt, epsxcmt, vxcmt)
    end

    if spinpol
        potxcir!(rhoir, magir, epsxcir, vxcir, bxcir)
    else
        potxcir!(rhoir, epsxcir, vxcir)
    end

    # Symmetrize
    symrfmt!(atoms, mt_vars, sym_vars, vxcmt)
    symrfir!(pw, sym_vars, vxcir)
    if spinpol
        symrfmt!(atoms, mt_vars, sym_vars, bxcmt)
        symrfir!(pw, sym_vars, bxcir)
    end

    # effective potential from sum of Coulomb and exchange-correlation potentials
    for ia in 1:Natoms
        @views vsmt[ia][:] .= vclmt[ia][:] .+ vxcmt[ia][:]
    end
    vsir[:] = vclir[:] + vxcir[:]

    # smoothing vsir is skipped (default is zero)
    
    # Generate the effective magnetic fields
    if spinpol
        genbs!( atoms, mt_vars, cfunir, ncmag,
            bfcmt, bfieldc, bxcmt, bxcir, bsmt, bsir
        )
    end

    # generating the tau-DFT effective potential is skipped

    return

end