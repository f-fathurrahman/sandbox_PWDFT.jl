function olpaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, O)
    
    atm2species = atoms.atm2species
    Ngwk = pw.gvecw.Ngw[ik]
    lmoapw = apwlo_vars.lmoapw
    apword = apwlo_vars.apword
    lmaxapw = mt_vars.lmaxapw

    isp = atm2species[ia]
    lmo = lmoapw[isp]
    a = zeros(ComplexF64, lmo, Ngwk)
    i = 0
    lm = 0
    for l in 0:lmaxapw, m in -l:l
        lm += 1
        for io in 1:apword[isp][l]
            i += 1
            a[i,1:Ngwk] .= apwalm[ia][1:Ngwk,io,lm]
        end
    end
    zmctmu!(a, a, O)
    return
end
