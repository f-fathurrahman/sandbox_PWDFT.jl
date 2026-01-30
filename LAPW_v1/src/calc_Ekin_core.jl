
function calc_Ekin_core(atoms, atsp_vars, core_states, mt_vars, vsmt)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    
    nstsp = atsp_vars.nstsp
    spcore = atsp_vars.spcore

    occcr = core_states.occcr
    evalcr = core_states.evalcr
    spincore = core_states.spincore
    rhocr = core_states.rhocr

    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    npmt = mt_vars.npmt
    npmtmax = maximum(mt_vars.npmt)
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo

    rfmt = zeros(Float64, npmtmax)

    # calculate the kinetic energy for core states
    Ekin_core = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        # sum of core eigenvalues
        for ist in 1:nstsp[isp]
            if spcore[isp][ist]
                Ekin_core += occcr[ia][ist] * evalcr[ia][ist]
            end
        end
        # core density
        rfmt[1:npmt[isp]] .= 0.0
        if spincore
            # spin-polarized core
            i = 1
            for ir in 1:nri
                rfmt[i] = rhocr[ir,ias,1] + rhocr[ir,ias,2]
                i += lmmaxi
            end
            for ir in (nri+1):nr
                rfmt[i] = rhocr[ir,ias,1] + rhocr[ir,ias,2]
                i += lmmaxo
            end
        else
            # spin-unpolarized core
            i = 1
            for ir in 1:nri
                rfmt[i] = rhocr[ia][ir,1]
                i += lmmaxi
            end
            for ir in (nri+1):nr
                rfmt[i] = rhocr[ia][ir,1]
                i += lmmaxo
            end
        end
        Ekin_core -= rf_mt_inner_prod(isp, mt_vars, rfmt[1:npmt[isp]], vsmt[ia])
    end
    return Ekin_core
end
