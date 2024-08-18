function hmllolo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, hlolo, H)

    atm2species = atoms.atm2species
    
    Ngwk = pw.gvecw.Ngw[ik]
    
    lmaxo = mt_vars.lmaxo
    idxlm = mt_vars.idxlm
    gntyry = mt_vars.gntyry

    nlorb = apwlo_vars.nlorb
    lorbl = apwlo_vars.lorbl
    idxlo = apwlo_vars.idxlo

    isp = atm2species[ia]

    for jlo in 1:nlorb[isp]
        l3 = lorbl[isp][jlo]
        for m3 in -l3:l3
            lm3 = idxlm[l3,m3]
            j = Ngwk + idxlo[ia][jlo][lm3]
            for ilo in 1:nlorb[isp]
                l1 = lorbl[isp][ilo]
                for m1 in -l1:l1
                    lm1 = idxlm[l1,m1]
                    i = Ngwk + idxlo[ia][ilo][lm1]
                    if i <= j
                        z1 = 0.0 + im*0.0
                        for l2 in 0:lmaxo
                            if mod(l1+l2+l3, 2) == 0
                                for m2 in -l2:l2
                                    lm2 = idxlm[l2,m2]
                                    z1 += gntyry[lm2,lm3,lm1] * hlolo[ia][lm2,jlo,ilo]
                                end
                            end
                        end
                        H[i,j] += z1
                    end # if
                end # m1
            end
        end
    end
    return 
end
