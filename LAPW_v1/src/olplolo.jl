function olplolo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, ololo, O)
    atm2species = atoms.atm2species
    Ngwk = pw.gvecw.Ngw[ik]
    idxlm = mt_vars.idxlm
    nlorb = apwlo_vars.nlorb
    idxlo = apwlo_vars.idxlo
    lorbl = apwlo_vars.lorbl
    #
    isp = atm2species[ia]
    for ilo in 1:nlorb[isp]
        l = lorbl[isp][ilo]
        for jlo in 1:nlorb[isp]
            if lorbl[isp][jlo] == l 
                for m in -l:l
                    lm = idxlm[l,m]
                    i = Ngwk + idxlo[ia][ilo][lm]
                    j = Ngwk + idxlo[ia][jlo][lm]
                    if i <= j 
                        O[i,j] += ololo[ia][ilo,jlo]
                    end # if 
                end # m
            end # if
        end # jlo
    end # ilo
    return
end

