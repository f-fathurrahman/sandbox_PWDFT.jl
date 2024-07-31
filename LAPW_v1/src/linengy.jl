function linengy!()
    #=
  USE m_atoms, ONLY: natmmax, idxas, natoms, nspecies
  USE m_apwlo, ONLY: nlorb, lorbord, lorbe, apword, apwe, &
               lorbl, lorbe0, lorbe, lorbve, apwe0, apwve, dlefe, &
               autolinengy, demaxbnd, epsband
  USE m_symmetry, ONLY: eqatoms
  USE m_muffin_tins, ONLY: nrmtmax, nrmti, nrmt, lmmaxi, lmmaxo, lmaxapw, rlmt
  USE m_constants, ONLY: y00, solsc
  USE m_density_pot_xc, ONLY: vsmt
  USE m_convergence, ONLY: iscl
  USE m_states, ONLY: efermi
    =#

    nnf = 0
    done = zeros(Bool, Natoms) # falses
    #
    nrmtmax = maximum(nrmt)
    vr = zeros(Float64, nrmtmax)

    # begin loops over atoms and species
    for ia in 1:Natoms
        #
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        iro = nri + 1
        #
        if done[ia]
            continue
        end
        i = 1
        for ir in 1:nri
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxi
        end
        for ir in iro:nr
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxo
        end 
        # APW functions
        for l in 0:lmaxapw
            for io in 1:apword[isp][l]
                if apwve[isp][l][io]
                    # check if previous radial functions have same default energies
                    for jo in 1:(io-1)
                        if apwve[isp][l][jo] 
                            if abs(apwe0[isp][l][io] - apwe0[isp][l][jo]) < 1.e-4
                                apwe[ia][l][io] = apwe[ia][l][jo]
                                @goto LABEL10 # next order io
                            end
                        end
                    end
                    # find the band energy starting from default
                    apwe[ia][l][io] = apwe0[isp][l][io]
                    @views rr = rlmt[isp][:,1]
                    fnd, apwe[ia][l][io] = findband!(l, rr, vr, apwe[ia][l][io])
                    if !fnd
                        nnf += 1
                    else 
                        # set linearization energy automatically if allowed/enabled
                        if autolinengy
                            apwe[ia][l][io] = efermi + dlefe
                        end
                    end
                end
                @LABEL10
            end # io order
        end
        #
        # local-orbital functions
        #
        for ilo in 1:nlorb[isp]
            for io in 1:lorbord[isp][ilo]
                if lorbve[isp][ilo][io] 
                    # check if previous radial functions have same default energies
                    for jo in 1:(io-1)
                        if lorbve[isp][ilo][jo]
                            if abs(lorbe0[isp][ilo][io] - lorbe0[isp][ilo][jo]) < 1.e-4
                                lorbe[ia][ilo][io] = lorbe[ia][ilo][jo]
                                @goto LABEL20 # next order
                            end
                        end
                    end
                    l = lorbl[isp][ilo]
                    # find the band energy starting from default
                    lorbe[ia][ilo][io] = lorbe0[isp][ilo][io]
                    @views rr = rlmt[isp][:,1]
                    fnd, lorbe[ia][ilo][io] = findband!(l, rr, vr, lorbe[ia][ilo][io])
                    if !fnd
                        nnf += 1
                    end
                else 
                    # set linearization energy automatically if allowed/enabled
                    if autolinengy
                        lorbe[ia][ilo][io] = efermi + dlefe
                    end
                end 
                @LABEL20 # next order
            end # over all order
        end # over all nlorb
        #
        done[ia] = true
        # copy to equivalent atoms
        for ja in 1:Natoms
            if !done[ja] && eqatoms[ia,ja]
                for l in 0:lmaxapw
                    for io in 1:apword[isp][l]
                        apwe[ja][l][io] = apwe[ia][l][io]
                    end 
                end
                for ilo in 1:nlorb[isp]
                    for io in 1:lorbord[isp][ilo]
                        lorbe[ja][ilo][io] = lorbe[ia][ilo][io]
                    end 
                end
                done[ja] = true
            end 
        end
    end # end loops over atoms
    #
    if nnf > 0
        @warn "Warning(linengy): could not find $nnf linearization energies"
    end
end

