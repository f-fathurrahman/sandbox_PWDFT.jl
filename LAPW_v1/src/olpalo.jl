function olpalo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, oalo, O)
    atm2species = atoms.atm2species
    Ngwk = pw.gvecw.Ngw[ik]
    idxlm = mt_vars.idxlm
    idxlo = apwlo_vars.idxlo
    lorbl = apwlo_vars.lorbl
    nlorb = apwlo_vars.nlorb
    apword = apwlo_vars.apword

    nmatk = size(O, 1)
    isp = atm2species[ia]
    for ilo in 1:nlorb[isp]
        l = lorbl[isp][ilo]
        for m in -l:l
            lm = idxlm[l,m]
            j = Ngwk + idxlo[ia][ilo][lm]
            k = (j-1)*nmatk
            for i in 1:Ngwk
                k += 1
                for io in 1:apword[isp][l]
                    O[k] += conj(apwalm[ia][i,io,lm]) * oalo[ia][io,ilo]
                end
            end
        end
    end
    return
end
