function rhocore!(atoms, mt_vars, elec_chgst, core_states, rhomt; magmt=nothing)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    nrmtmax = maximum(nrmt)
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    wrmt = mt_vars.wrmt

    fr = zeros(Float64, nrmtmax)
    y00 = 0.28209479177387814347

    rhocr = core_states.rhocr
    spincore = core_states.spincore
    nspncr = core_states.nspncr
    
    chgcr = elec_chgst.chgcr
    chgcrlk = elec_chgst.chgcrlk

    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        iro = nri + 1
        ss = 0.0
        # loop over spin channels
        for ispn in 1:nspncr
            # add the core density to the muffin-tin density
            i = 1
            for ir in 1:nri
                rhomt[ia][i] += rhocr[ia][ir,ispn]
                i += lmmaxi
            end
            for ir in iro:nr
                rhomt[ia][i] += rhocr[ia][ir,ispn]
                i += lmmaxo
            end
            # compute the core charge inside the muffin-tins
            t1 = dot(wrmt[isp], rhocr[ia][:,ispn])
            ss += 4π * y00 * t1
        end
        # core leakage charge
        chgcrlk[ia] = chgcr[isp] - ss
        # add to the magnetization in the case of a spin-polarized core
        if spincore
            # compute the moment in the muffin-tin
            for idm in 1:ndmag
                rfmtlm(1,nr,nri,magmt(:,ias,idm),fr)
                rf_mt_lm!(1, isp, mt_vars, magmt[ia][:,idm], fr)
                t1 = dot(wrmt[isp], fr[1:nr])
                v[idm] = 4π * y00 * t1
            end
            # normalize
            if ncmag
                t1 = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
            else
                t1 = abs(v(1))
            end
            if t1 > 1e-10
                v[:] *= 1/t1
            end
            # add the core magnetization to the total
            i = 1
            for ir in 1:nri
                t1 = abs( rhocr[ia][ir,1] - rhocr[ia][ir,2] )
                @. magmt[ia][i,:] += t1*v[:]
                i += lmmaxi
            end
            for ir in iro:nr
                t1 = abs( rhocr[ia][ir,1] - rhocr[ia][ir,2] )
                @. magmt[ia][i,:] += t1*v[:]
                i += lmmaxo
            end
        end # spincore
    end
    return
end
