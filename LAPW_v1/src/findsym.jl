function findsym!(atoms, apl1, apl2, nsym, lspl, lspn, iea; epslat=1e-6)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    jea = zeros(Int64, Natoms)
    sl = zeros(Float64, 3, 3)

    nsym = 0
    sl = zeros(Float64, 3, 3)
    # loop over lattice symmetries (spatial rotations)
    for isym in 1:nsymlat
        # make real copy of lattice rotation symmetry
        @views sl[:,:] = Float64.(symlat[isym][:,:])
        # loop over species
        for ia in 1:Natoms
            isp = atm2species[ia]
            # map apl1 coordinates to [0,1) and store in apl3
            apl3[:,ia] = apl1[:,ia]
            @views r3frac!(apl3[:,ia], epslat=epslat)
            #
            for ja in 1:Natoms
                if atm2species[ja] != isp
                    continue
                end
                # apply lattice symmetry to atomic positions
                @views v[:] = sl[:,1]*apl2[1,ja] + sl[:,2]*apl2[2,ja] + sl[:,3]*apl2[3,ja]
                # map coordinates to [0,1)
                r3frac!(v, epslat=epslat)
                # check if atomic positions are invariant
                t1 = abs(apl3[1,ia]-v[1]) + abs(apl3[2,ia]-v[2]) + abs(apl3[3,ia]-v[3])
                if t1 < epslat 
                    # equivalent atom index
                    jea[ia] = ja
                    @goto LABEL10 # continue ?
                end
                # not invariant so try new spatial rotation
                @goto LABEL40
                @label LABEL10 #10 CONTINUE
            end
        end
    
        # all atomic positions invariant at this point
        jsym = 1

        # 
        # .... spin polarized stuffs is removed
        # spinpol    
    
        # everything invariant so add symmetry to set
        nsym += 1
        lspl[nsym] = isym
        lspn[nsym] = jsym
        for ia in 1:Natoms
            iea[ia,nsym] = jea[ia]
        end

        @label LABEL40 # continue
    
    end # ! end loop over spatial rotations 

    return
end

