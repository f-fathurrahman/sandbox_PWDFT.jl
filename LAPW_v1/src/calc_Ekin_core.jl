
function calc_Ekin_core()

    # calculate the kinetic energy for core states
    Ekin_core = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        # sum of core eigenvalues
        for ist in 1:nstsp[isp]
            if spcore[isp][ist]
                Ekin_core += occcr[ist,ia] * evalcr[ist,ia]
            end
        end
        # core density
        rfmt[1:npmt[isp]] = 0.0
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
                rfmt[i] = rhocr[ir,ias,1]
                i += lmmaxi
            end
            for ir in (nri+1):nr
                rfmt[i] = rhocr[ir,ias,1]
                i += lmmaxo
            end
        end
        Ekin_core -= rf_mt_inner_prod(nr, nri, wrmt[isp], rfmt, vsmt[ia])
    end
    return Ekin_core
end
