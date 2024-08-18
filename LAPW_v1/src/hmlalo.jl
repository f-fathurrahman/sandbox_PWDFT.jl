function hmlalo!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, hloa, H)

    atm2species = atoms.atm2species
    
    Ngwk = pw.gvecw.Ngw[ik]
    
    lmaxapw = mt_vars.lmaxapw
    lmaxo = mt_vars.lmaxo
    idxlm = mt_vars.idxlm
    gntyry = mt_vars.gntyry

    apword = apwlo_vars.apword
    nlorb = apwlo_vars.nlorb
    lorbl = apwlo_vars.lorbl
    idxlo = apwlo_vars.idxlo

    nmatk = size(H, 1)
    isp = atm2species[ia]

    for ilo in 1:nlorb[isp]
        l1 = lorbl[isp][ilo]
        for m1 in -l1:l1
            lm1 = idxlm[l1,m1]
            j = Ngwk + idxlo[ia][ilo][lm1]
            lm3 = 0
            for l3 in 0:lmaxapw, m3 in -l3:l3
                lm3 = lm3 + 1
                for io in 1:apword[isp][l3]
                    z1 = 0.0 + im*0.0
                    for l2 in 0:lmaxo
                        if mod(l1+l2+l3, 2) == 0
                            for m2 in -l2:l2
                                lm2 = idxlm[l2,m2]
                                z1 += gntyry[lm2,lm3,lm1] * hloa[ia][lm2,io,l3,ilo]
                            end
                        end # if
                    end 
                    # note that what is actually computed is the Hermitian conjugate of <lo|H|APW>
                    if abs(real(z1)) + abs(imag(z1)) > 1.e-14
                        k = (j-1)*nmatk
                        for i in 1:Ngwk
                            k += 1
                            H[k] += conj(z1 * apwalm[ia][i,io,lm3] )
                        end
                    end # if
                end # io
            end 
        end
    end
    return

end
