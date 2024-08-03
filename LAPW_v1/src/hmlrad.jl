function hmlrad()

    fr = zeros(Float64, nrmtmax)
  
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        iro = nri + 1
        npi = npmti[isp]
        #---------------------------
        #     APW-APW integrals     
        #---------------------------
        for l1 in 0:lmaxapw
            for io in 1:apword[isp][l1]
                for l3 in 0:lmaxapw
                    for jo in 1:apword[isp][l3]
                        if l1 == l3
                            fr[1:nr] .= apwfr[ia][l1][io][1:nr,1] .* apwfr[ia][l3][jo][1:nr,2]
                            t1 = dot( wrmt[isp][1:nr], fr[1:nr] )
                            haa[ia][1,jo,l3,io,l1] = t1/y00
                        else 
                            haa[ia][1,jo,l3,io,l1] = 0.0
                        end 
                        #
                        if l1 >= l3
                            #
                            lm2 = 1
                            #
                            for l2 in 1:lmaxi
                                for m2 in -l2:l2
                                    lm2 = lm2 + 1
                                    i = lm2
                                    for ir in 1:nri
                                        t1 = apwfr[ia][l1][io][ir,1] * apwfr[ia][l3][jo][ir,1]
                                        fr[ir] = t1*vsmt[ia][i]
                                        i = i + lmmaxi
                                    end 
                                    for ir in iro:nr
                                        t1 = apwfr[ia][l1][io][ir,1] * apwfr[ia][l3][jo][ir,1]
                                        fr[ir] = t1*vsmt[ia][i]
                                        i = i + lmmaxo
                                    end 
                                    t1 = dot( wrmt[isp][1:nr], fr[1:nr] )
                                    haa[ia][lm2,jo,l3,io,l1] = t1
                                    haa[ia][lm2,io,l1,jo,l3] = t1
                                end 
                            end
                            # 
                            for l2 in (lmaxi+1):lmaxo
                                for m2 in -l2:l2
                                    lm2 = lm2 + 1
                                    i = npi + lm2
                                    for ir in iro:nr
                                        t1 = apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)
                                        fr[ir] = t1*vsmt[ia][i]
                                        i = i + lmmaxo
                                    end 
                                    t1 = dot( wrmt[isp][iro:nr], fr[iro:nr] )
                                    haa[ia][lm2,jo,l3,io,l1] = t1
                                    haa[ia][lm2,io,l1,jo,l3] = t1
                                end 
                            end 
                        end # if 
                    end 
                end 
            end 
        end 
        #-------------------------------------
        #     local-orbital-APW integrals     
        #-------------------------------------
        for ilo in 1:nlorb[isp]
            l1 = lorbl[isp][ilo]
            for l3 in 0:lmaxapw
                for io in 1:apword[isp][l3]
                    if l1 == l3
                        fr[1:nr] .= lofr[ia][ilo][1:nr,1] .* apwfr[ia][l3][io][1:nr,2]
                        t1 = dot( wrmt[isp][1:nr], fr[1:nr] )
                        hloa[ia][1,io,l3,ilo] = t1/y00
                    else 
                        hloa[ia][1,io,l3,ilo] = 0.0
                    end # if
                    lm2 = 1
                    for l2 in 1:lmaxi
                        for m2 in -l2:l2
                            lm2 = lm2 + 1
                            i = lm2
                            for ir in 1:nri
                                t1 = lofr[ia][ilo][ir,1] * apwfr[ia][l3][io][ir,1]
                                fr[ir] = t1*vsmt[ia][i]
                                i = i+lmmaxi
                            end
                            for ir in (nri+1):nr
                                t1 = lofr[ia][ilo][ir,1] * apwfr[ia][l3][io][ir,1]
                                fr[ir] = t1*vsmt[ia][i]
                                i = i + lmmaxo
                            end 
                            t1 = dot( wrmt[isp][1:nr], fr[1:nr] )
                            hloa[ia][lm2,io,l3,ilo] = t1
                        end 
                    end
                    for l2 in (lmaxi+1):lmaxo
                        for m2 in -l2:l2
                            lm2 = lm2 + 1
                            i = npi + lm2
                            for ir in iro:nr
                                t1 = lofr[ia][ilo][ir,1] * apwfr[ia][l3][ir,1]
                                fr[ir] = t1*vsmt[ia][i]
                                i = i + lmmaxo
                            end 
                            t1 = dot( wrmt[isp][iro:nr], fr[iro:nr] )
                            hloa[ia][lm2,io,l3,ilo] = t1
                        end 
                    end
                end 
            end
        end
        #-----------------------------------------------
        #     local-orbital-local-orbital integrals     
        #-----------------------------------------------
        for ilo in 1:nlorb[isp]
            l1 = lorbl[isp][ilo]
            for jlo in 1:nlorb[isp]
                l3 = lorbl[isp][jlo]
                if l1 == l3
                    fr[1:nr] = lofr[ia][ilo][1:nr,1] .* lofr[ia][jlo][1:nr,2]
                    t1 = dot( wrmt[isp][1:nr], fr[1:nr] )
                    hlolo[ia][1,jlo,ilo] = t1/y00
                else 
                    hlolo[ia][1,jlo,ilo] = 0.0
                end 
                lm2 = 1
                for l2 in 1:lmaxi
                    for m2 in -l2:l2
                        lm2 = lm2 + 1
                        i = lm2
                        for ir in 1:nri
                            t1 = lofr[ia][ilo][ir,1] * lofr[ia][jlo][ir,1]
                            fr[ir] = t1*vsmt[ia][i]
                            i = i + lmmaxi
                        end 
                        for ir in iro:nr
                            t1 = lofr[ia][ilo][ir,1] * lofr[ia][jlo][ir,1]
                            fr[ir] = t1*vsmt[ia][i]
                            i = i + lmmaxo
                        end 
                        t1 = dot( wrmt[isp][1:nr], fr[1:nr] )
                        hlolo[ia][lm2,jlo,ilo] = t1
                    end 
                end 
                for l2 in (lmaxi+1):lmaxo
                    for m2 in -l2:l2
                        lm2 = lm2 + 1
                        i = npi + lm2
                        for ir in iro:nr
                            t1 = lofr[ia][ilo][ir,1] * lofr[ia][jlo][ir,1]
                            fr[ir] = t1*vsmt[ia][i]
                            i = i + lmmaxo
                        end
                        t1 = dot( wrmt[isp][iro:nr], fr[iro:nr] )
                        hlolo[ia][lm2,jlo,ilo] = t1
                    end 
                end 
            end 
        end 
    end
    return
end