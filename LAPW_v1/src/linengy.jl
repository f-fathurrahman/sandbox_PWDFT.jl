function linengy!(atoms, eqatoms, mt_vars, vsmt, efermi::Float64, apwlo_vars)
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

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    lmmaxo = mt_vars.lmmaxo
    lmmaxi = mt_vars.lmmaxi
    lmaxapw = mt_vars.lmaxapw
    rlmt = mt_vars.rlmt

    apwe = apwlo_vars.apwe
    apword = apwlo_vars.apword
    apwve = apwlo_vars.apwve
    nlorb = apwlo_vars.nlorb
    lorbve = apwlo_vars.lorbve
    lorbord = apwlo_vars.lorbord
    autolinengy = apwlo_vars.autolinengy
    lorbl = apwlo_vars.lorbl
    lorbe0 = apwlo_vars.lorbe0
    lorbe = apwlo_vars.lorbe
    apwe0 = apwlo_vars.apwe0

    nnf = 0
    done = zeros(Bool, Natoms) # falses
    #
    nrmtmax = maximum(nrmt)
    vr = zeros(Float64, nrmtmax)

    y00 = 0.28209479177387814347

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
        @info "lmaxapw = $lmaxapw"
        for l in 0:lmaxapw
            for io in 1:apword[isp][l]
                if apwve[isp][l][io]
                    @info "Enter here 69"
                    # check if previous radial functions have same default energies
                    for jo in 1:(io-1)
                        if apwve[isp][l][jo] 
                            if abs(apwe0[isp][l][io] - apwe0[isp][l][jo]) < 1.e-4
                                @info "Same default energies for apw"
                                apwe[ia][l][io] = apwe[ia][l][jo]
                                @goto LABEL10 # next order io
                            end
                        end
                    end
                    # find the band energy starting from default
                    apwe[ia][l][io] = apwe0[isp][l][io]
                    @views rr = rlmt[isp][:,1]
                    @info "Calling findband for APW"
                    fnd, apwe[ia][l][io] = findband!(l, rr, vr, apwe[ia][l][io])
                    @info "fnd = $fnd"
                    if !fnd
                        nnf += 1
                    else 
                        # set linearization energy automatically if allowed/enabled
                        if autolinengy
                            apwe[ia][l][io] = efermi + dlefe
                        end
                    end
                end  # if apwve
                @label LABEL10
            end # io order
        end
        #
        # local-orbital functions
        #
        for ilo in 1:nlorb[isp]
            for io in 1:lorbord[isp][ilo]
                if lorbve[isp][ilo][io]
                    println("lorbve = ", lorbve[isp][ilo][io])
                    @info "Enter here 104"
                    # check if previous radial functions have same default energies
                    for jo in 1:(io-1)
                        if lorbve[isp][ilo][jo]
                            if abs(lorbe0[isp][ilo][io] - lorbe0[isp][ilo][jo]) < 1.e-4
                                @info "Same default energies for local orbitals"
                                lorbe[ia][ilo][io] = lorbe[ia][ilo][jo]
                                @goto LABEL20 # next order
                            end
                        end
                    end
                    l = lorbl[isp][ilo]
                    # find the band energy starting from default
                    lorbe[ia][ilo][io] = lorbe0[isp][ilo][io]
                    @views rr = rlmt[isp][:,1]
                    @info "Calling findband for local orbitals"
                    fnd, lorbe[ia][ilo][io] = findband!(l, rr, vr, lorbe[ia][ilo][io])
                    @info "fnd = $(fnd)"
                    if !fnd
                        nnf += 1
                    end
                else 
                    # set linearization energy automatically if allowed/enabled
                    if autolinengy
                        lorbe[ia][ilo][io] = efermi + dlefe
                    end
                end  # if lorbve
                @label LABEL20 # next order
            end # over all order
        end # over all nlorb
        #
        done[ia] = true
        # copy to equivalent atoms
        for ja in 1:Natoms
            if !done[ja] && eqatoms[ia,ja]
                @info "Do for eqatoms"
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
    @info "nnf = $nnf"
    #
    if nnf > 0
        @warn "Warning(linengy): could not find $nnf linearization energies"
    end
end

