function gencore!(atoms, eqatoms, atsp_vars, mt_vars, vsmt, core_states)

    NspinCore = core_states.nspncr
    @assert NspinCore == 1

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    nrsp = atsp_vars.nrsp
    nstsp = atsp_vars.nstsp
    nrspmax = maximum(nrsp)
    nstspmax = maximum(nstsp)
    vrsp = atsp_vars.vrsp
    spcore = atsp_vars.spcore
    nsp = atsp_vars.nsp
    lsp = atsp_vars.lsp
    ksp = atsp_vars.ksp
    rsp = atsp_vars.rsp

    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    rlmt = mt_vars.rlmt

    done = zeros(Bool, Natoms)
    vr = zeros(Float64, nrspmax)
    evals = zeros(Float64, nstspmax)

    rhocr = core_states.rhocr
    evalcr = core_states.evalcr
    rwfcr = core_states.rwfcr
    occcr = core_states.occcr

    y00 = 1/(2*sqrt(pi)) # spherical harmonics Y_00

    # loop over all atoms
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        nrs = nrsp[isp]
        if done[ia]
            continue
        end
        #
        # ... Spin-polarized core is skipped ...
        #
        #
        # loop over spin channels
        for ispincr in 1:NspinCore
            # use the spherical part of the crystal Kohn-Sham potential
            @views rf_mt_lm!(1, isp, mt_vars, vsmt[ia], vr[1:nr])
            @views vr[1:nr] = vr[1:nr]*y00
            # spin-up and -down potentials for polarized core
            # ... spincore case is skipped ...
            #
            # append the Kohn-Sham potential from the atomic calculation for r > R_MT
            t1 = vr[nr] - vrsp[isp][nr]
            for ir in (nr+1):nrs
                vr[ir] = vrsp[isp][ir] + t1
            end
            # ispin_core index is not implemented
            rhocr[ia][1:nr,ispincr] .= 0.0
            for ist in 1:nstsp[isp]
                # skip if this state is not core state
                if !spcore[isp][ist]
                    continue
                end
                # solve the Dirac equation
                evals[ist] = evalcr[ia][ist]
                @views psi1 = rwfcr[ia][:,1,ist]
                @views psi2 = rwfcr[ia][:,2,ist]
                𝓃 = nsp[isp][ist]
                𝓁 = lsp[isp][ist]
                𝓀 = ksp[isp][ist]
                #@info "size rsp[isp] = $(size(rsp[isp]))"
                #@info "size vr = $(size(vr))"
                evals[ist] = rdirac!( 𝓃, 𝓁, 𝓀, rsp[isp], vr, evals[ist], psi1, psi2 )
                #
                # ... spincore case is skipped ...
                #
                evalcr[ia][ist] = evals[ist]
                t1 = occcr[ia][ist]
                # add to the core density
                for ir in 1:nr
                    rhocr[ia][ir,ispincr] += t1*( rwfcr[ia][ir,1,ist]^2 + rwfcr[ia][ir,2,ist]^2 )
                end
            end
            for ir in 1:nr
                rhocr[ia][ir,ispincr] *= ( rlmt[isp][ir,-2]*y00 )
            end
        end # end loop over spin channels
        #
        done[ia] = true
        #
        # copy to equivalent atoms
        for ja in 1:Natoms
            if !done[ja] && eqatoms[ia,ja]
                for ist in 1:nstsp[isp]
                    # skip if not a core state
                    if !spcore[isp][ist]
                        continue
                    end
                    evalcr[ja][ist] = evalcr[ia][ist]
                    rwfcr[ja][1:nrs,:,ist] .= rwfcr[ia][1:nrs,:,ist]
                end
                for ispincr in 1:NspinCore
                    rhocr[ja][1:nr,ispincr] .= rhocr[ia][1:nr,ispincr]
                end
                done[ja] = true
            end
        end
    end # end loop over species and atoms

    return
end