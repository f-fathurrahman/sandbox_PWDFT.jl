function occupy!(
    kpoints, apwlo_vars, elec_chgst;
    NiterMax=1000, epsocc=1e-8
)

    evalsv = elec_chgst.evalsv
    swidth = elec_chgst.swidth
    efermi = elec_chgst.efermi
    occsv = elec_chgst.occsv
    chgval = elec_chgst.chgval

    if elec_chgst.nspinor == 2
        occmax = 1.0
    elseif elec_chgst.nspinor == 1
        occmax = 2.0
    else
        @error("Invalid nspinor value")
    end

    Nkpt = kpoints.Nkpt
    wkpt = kpoints.wk
    e0min = apwlo_vars.e0min

    nstsv = elec_chgst.nstsv

    # FIXME: Is this enough?
    @assert Nkpt == size(evalsv, 2)

    # XXX: automatic smearing width is disabled
    
    # find minimum and maximum eigenvalues
    e0 = evalsv[1,1]
    e1 = e0
    for ik in 1:Nkpt
        for ist in 1:nstsv
            e = evalsv[ist,ik]
            if e < e0
                e0 = e
            end
            if e > e1
                e1 = e
            end
        end
    end
    
    if e0 < e0min
        println("WARNING: minimum eigenvalue less than minimum linearization energy : ", e0, e0min)
    end
    
    t1 = 1.0/swidth
    # determine the Fermi energy using the bisection method
    iterOcc = 0
    while iterOcc <= NiterMax
        elec_chgst.efermi = 0.50*(e0 + e1)
        chg = 0.0
        for ik in 1:Nkpt
            for ist in 1:nstsv
                e = evalsv[ist,ik]
                if e < e0min
                    occsv[ist,ik] = 0.0
                else
                    x = (elec_chgst.efermi - e)*t1
                    occsv[ist,ik] = occmax*stheta_fd(x) # smearing
                    chg += wkpt[ik]*occsv[ist,ik]
                end 
            end
        end
        if chg < chgval
            e0 = elec_chgst.efermi
        else
            e1 = elec_chgst.efermi
        end
        if abs(e1-e0) < 1.e-12
            break
        end
        iterOcc += 1
    end
    if iterOcc > NiterMax 
        println("WARNING: could not find Fermi energy")
    end


    # find the density of states at the Fermi surface in units of
    # states/Hartree/unit cell
    fermidos = 0.0
    for ik in 1:Nkpt
        for ist in 1:nstsv
            x = (evalsv[ist,ik] - elec_chgst.efermi)*t1
            fermidos += wkpt[ik]*sdelta_fd(x)*t1
        end 
        if abs(occsv[nstsv,ik]) > epsocc
            println("WARNING: not enough empty states for k-point ", ik)
        end
    end
    fermidos = fermidos*occmax
    #XXX Need this?

    return
end
