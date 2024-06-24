function findsym!(
    sym_vars::SymmetryVars,
    atoms::Atoms,
    apl1, apl2, lspl, lspn, iea;
    epslat=1e-6
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nsymlat = sym_vars.nsymlat
    symlat = sym_vars.symlat

    jea = zeros(Int64, Natoms)
    sl = zeros(Float64, 3, 3)
    apl3 = zeros(Float64, 3, Natoms)
    v = zeros(Float64, 3)

    nsym = 0
    sl = zeros(Float64, 3, 3)
    # loop over lattice symmetries (spatial rotations)
    for isym in 1:nsymlat
        #@info "Try isym = $isym"

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
                print("Loop ia=$ia ja=$ja isym=$isym")
                if atm2species[ja] != isp
                    @info "Skipping this atom index"
                    continue
                end
                # apply lattice symmetry to atomic positions
                @views v[:] = sl[:,1]*apl2[1,ja] + sl[:,2]*apl2[2,ja] + sl[:,3]*apl2[3,ja]
                # map coordinates to [0,1)
                r3frac!(v, epslat=epslat)
                # check if atomic positions are invariant
                t1 = abs(apl3[1,ia]-v[1]) + abs(apl3[2,ia]-v[2]) + abs(apl3[3,ia]-v[3])
                #@info "t1 = $t1"
                if t1 < epslat 
                    println("*** Equivalent atoms: ia=$ia ja=$ja isym=$isym")
                    # equivalent atom index
                    jea[ia] = ja
                    @goto LABEL10 # continue ?
                end
                # not invariant so try new spatial rotation
                println(" - Not invariant, trying new spatial rotation")
                @goto LABEL40 # jump to the end of loop over symlat
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
        println("nsym = ", nsym)

        @info "Trying new spatial rotation"
        @label LABEL40 # continue
    
    end # ! end loop over spatial rotations 

    @info "nsym = $nsym"
    return nsym

    return
end

