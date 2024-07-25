function symrfmt!(atoms, mt_vars, sym_vars, rfmt)
  
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    nsymcrys = sym_vars.nsymcrys
    symlatc = sym_vars.symlatc
    ieqatom = sym_vars.ieqatom
    lsplsymc = sym_vars.lsplsymc
    isymlat = sym_vars.isymlat

    t0 = 1.0/nsymcrys

    npmt = mt_vars.npmt
    npmtmax = maximum(npmt)
    rfmt1 = zeros(Float64, npmtmax, Natoms) # This is larger than necessary
    rfmt2 = zeros(Float64, npmtmax)

    done = zeros(Bool, Natoms)
    for isp in 1:Nspecies
        # make a copy of the input function
        for ia in 1:Natoms
            if isp != atm2species[ia]
                continue
            end
            @views rfmt1[1:npmt[isp],ia] .= rfmt[ia][:]
        end
        done[:] .= false
        # loop over atoms
        for ia in 1:Natoms
            #@info "Begin isp=$isp ia=$ia"
            #
            if isp != atm2species[ia]
                continue
            end
            #
            if done[ia]
                #@info "This ia=$ia is done, skipping"
                continue
            end
            # zero out the output
            rfmt[ia][1:npmt[isp]] .= 0.0
            #println("CHECK1: sum abs rfmt[ia] = ", sum(abs.(rfmt[ia])))
            # loop over crystal symmetries
            for isym in 1:nsymcrys
                #
                # index to spatial rotation lattice symmetry
                lspl = lsplsymc[isym]
                #
                # equivalent atom index (symmetry rotates atom ja into atom ia)
                ja = ieqatom[ia,isym]
                #
                # apply the rotation to the muffin-tin function
                @views rotrfmt!(mt_vars, isp, symlatc[lspl], rfmt1[:,ja], rfmt2)
                # XXX rfmt1 is a 2d array
                #
                # accumulate in original function array
                @views rfmt[ia][1:npmt[isp]] .+= rfmt2[1:npmt[isp]]
            end 
            # normalize by 1/nsymcrys
            rfmt[ia][1:npmt[isp]] .*= t0
            done[ia] = true # This atom has been processed
            #println("CHECK2: sum abs rfmt[ia] = ", sum(abs.(rfmt[ia])))
            # rotate into equivalent atoms
            for isym in 1:nsymcrys
                ja = ieqatom[ia,isym]
                if done[ja]
                    #println("This ja=$ja is done, skipping")
                    continue
                end
                #
                # inverse symmetry (which rotates atom ia into atom ja)
                lspl = isymlat[lsplsymc[isym]]
                #
                # rotate symmetrized function into equivalent muffin-tin
                rotrfmt!(mt_vars, isp, symlatc[lspl], rfmt[ia], rfmt[ja])
                #println("CHECK3: sum abs rfmt[ja] = ", sum(abs.(rfmt[ja])))
                done[ja] = true
            end
        end 
    end
    return
end

