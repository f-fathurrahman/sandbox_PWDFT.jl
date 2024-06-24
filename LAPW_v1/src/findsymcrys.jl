# will modify some fields of atoms
function findsymcrys!(
    sym_vars::SymmetryVars,
    atoms;
    tshift=true,
    epslat=1e-6,
    symtype=1
)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    LatVecs = atoms.LatVecs
    atm2species = atoms.atm2species
    atposc = atoms.positions
    atposl = inv(atoms.LatVecs)*atposc

    vtlsymc = sym_vars.vtlsymc
    lsplsymc = sym_vars.lsplsymc
    lspnsymc = sym_vars.lspnsymc
    vtcsymc = sym_vars.vtcsymc
    tv0symc = sym_vars.tv0symc

    symlat = sym_vars.symlat
    isymlat = sym_vars.isymlat

    MAX_SYM_CRYS = 172

    ieqatom = zeros(Int64, Natoms, MAX_SYM_CRYS)    
    eqatoms = zeros(Bool, Natoms, Natoms)

    # store position of first atom
    v0 = zeros(Float64, 3)
    if Natoms > 0
        @views v0[:] = atposl[:,1]
    end

    # find the smallest set of atoms
    # XXX: Need this?
    Natoms_per_species = zeros(Int64, Nspecies)
    for isp in 1:Nspecies
        Natoms_per_species[isp] = count(atm2species .== isp)
    end
    isp_smallest = 1
    for isp in 1:Nspecies
        if Natoms_per_species[isp] < Natoms_per_species[isp_smallest]
            isp_smallest = isp
        end
    end

    v1 = zeros(Float64, 3)
    v2 = zeros(Float64, 3)

    # atoms.positions will be shifted
    if tshift
        println("Unshifted atomic positions (before shifting)")
        for ia in 1:Natoms
            @printf("%18.10f %18.10f %18.10f\n", atposl[1,ia], atposl[2,ia], atposl[3,ia])
        end
        # shift basis so that the first atom in the smallest atom set is at the origin
        #
        # find the first atom with the species of smallest number of atoms
        for ia in 1:Natoms
            isp = atm2species[ia]
            if isp == isp_smallest
                v1[:] = atposl[:,ia]
                break
            end
        end
        println("v1 = ", v1)
        #
        for ia in 1:Natoms
            # shift atom
            atposl[:,ia] = atposl[:,ia] - v1[:]
            # map lattice coordinates back to [0,1)
            r3frac!(atposl[:,ia])
            # determine the new Cartesian coordinates
            atposc[:,ia] = LatVecs * atposl[:,ia]
        end 
        println("Shifted atomic positions")
        for ia in 1:Natoms
            @printf("%18.10f %18.10f %18.10f\n", atposl[1,ia], atposl[2,ia], atposl[3,ia])
        end
    end


    # determine possible translation vectors from smallest set of atoms

    #n = max( natoms(is)*natoms(is), 1 )
    @info "Natoms_per_species = $Natoms_per_species"
    @info "isp_smallest = $isp_smallest"

    vtl = zeros(Float64, 3, Natoms_per_species[isp_smallest]*Natoms_per_species[isp_smallest])
    n = 1
    vtl[:,1] .= 0.0 # not really needed
    for ia in 1:Natoms_per_species[isp_smallest]
        for ja in 2:Natoms_per_species[isp_smallest]
            # compute difference between two atom vectors
            v1[:] = atposl[:,ia] - atposl[:,ja]
            # map lattice coordinates to [0,1)
            r3frac!(v1)
            # 
            # .... tefield specific treatment is removed
            # 
            for i in 1:n
                t1 = abs( vtl[1,i] - v1[1] ) + abs( vtl[2,i] - v1[2] ) + abs( vtl[3,i] - v1[3] )
                if t1 < epslat
                    @goto LABEL10
                end
            end 
            n += 1
            vtl[:,n] = v1[:]
            @label LABEL10
        end
    end

    println()
    println("Translation vectors:")
    println("--------------------")
    for i in 1:n
        @printf("%18.10f %18.10f %18.10f\n", vtl[1,i], vtl[2,i], vtl[3,i])
    end

    # no translations required when symtype=0,2 (F. Cricchio)
    if symtype != 1
        @info "No translations are required for this symtype=$symtype"
        n = 1
    end

    println("n = ", n)  # XXX: change to Ntranslations

    eqatoms = zeros(Bool, Natoms, Natoms)
    apl = zeros(Float64, 3, Natoms)
    nsymcrys = 0
    lspl = zeros(Int64, 48)
    lspn = zeros(Int64, 48)
    iea = zeros(Int64, Natoms, 48)

    # loop over all possible translations
    for i in 1:n
        # construct new array with translated positions
        for ia in 1:Natoms
            apl[:,ia] = atposl[:,ia] + vtl[:,i]
        end
        # find the symmetries for current translation
        nsym = findsym!(sym_vars, atoms, atposl, apl, lspl, lspn, iea)
        #
        for isym in 1:nsym
            nsymcrys = nsymcrys + 1
            if nsymcrys > MAX_SYM_CRYS 
                error("Too many nsymcrys = $(nsymcrys)")
            end
            @views vtlsymc[:,nsymcrys] = vtl[:,i]
            lsplsymc[nsymcrys] = lspl[isym]
            lspnsymc[nsymcrys] = lspn[isym]
            for ia in 1:Natoms
                ja = iea[ia,isym]
                ieqatom[ia,nsymcrys] = ja
                eqatoms[ia,ja] = true
                eqatoms[ja,ia] = true
            end 
        end
    end

    @info "nsymcrys = $nsymcrys"

    tsyminv = false
    for isym in 1:nsymcrys
        i = lsplsymc[isym]
        # check if inversion symmetry is present
        if symlat[i] == -symlat[1]
            tsyminv = true
            # make inversion the second symmetry element (the identity is the first)
            v1[:] = vtlsymc[:,isym]; vtlsymc[:,isym] = vtlsymc[:,2]; vtlsymc[:,2] = v1[:]
            i = lsplsymc[isym]; lsplsymc[isym] = lsplsymc[2]; lsplsymc[2] = i
            i = lspnsymc[isym]; lspnsymc[isym] = lspnsymc[2]; lspnsymc[2] = i
            for ia in 1:Natoms
                i = ieqatom[ia,isym]
                ieqatom[ia,isym] = ieqatom[ia,2]
                ieqatom[ia,2] = i
            end
            break
        end
    end


    # if inversion exists THEN  shift basis so that inversion center is at origin
    if tsyminv && tshift 
        v1[:] = v1[:]/2.0
        for ia in 1:Natoms
            isp = atm2species[ia] # XXX not used
            # shift atom
            atposl[:,ia] = atposl[:,ia] + v1[:]
            # map lattice coordinates back to [0,1)
            @views r3frac!(atposl[:,ia], epslat=epslat)
            # map lattice coordinates to [-0.5,0.5)
            for i in 1:3
                if atposl[i,ia] > 0.5
                    atposl[i,ia] = atposl[i,ia] - 1.0
                end
            end
            # determine the new Cartesian coordinates
            @views r3mv!(LatVecs, atposl[:,ia], atposc[:,ia])
            # @views atposc[:,ia] = avec * atposl[:,ia]
        end # ia
        #
        # recalculate crystal symmetry translation vectors
        for isym in 1:nsymcrys
            ilspl = isymlat[lsplsymc[isym]]
            ss = symlat[ilspl] # shorthand
            v2[:] = ss[:,1]*v1[1] + ss[:,2]*v1[2] + ss[:,3]*v1[3]
            vtlsymc[:,isym] = vtlsymc[:,isym] - v1[:] + v2[:]
            @views r3frac!(vtlsymc[:,isym], epslat=epslat)
        end
    end # if

    # translation vector in Cartesian coordinates
    for isym in 1:nsymcrys
        @views r3mv!(LatVecs, vtlsymc[:,isym], vtcsymc[:,isym])
    end

    # set flag for zero translation vector
    for isym in 1:nsymcrys
        t1 = abs(vtlsymc[1,isym]) + abs(vtlsymc[2,isym]) + abs(vtlsymc[3,isym])
        if t1 < epslat 
            tv0symc[isym] = true
        else
            tv0symc[isym] = false
        end 
    end

    # check inversion does not include a translation
    if tsyminv 
        if !tv0symc[2]
            tsyminv = false
        end
    end 

    if Natoms > 0 
        v1[:] = atposl[:,1] - v0[:]
        t1 = abs(v1[1]) + abs(v1[2]) + abs(v1[3])
        if t1 > epslat
            println("INFO: atomic basis shift (latttice)")
            println(v1)
        end
    end

    sym_vars.ieqatom = ieqatom
    sym_vars.eqatoms = eqatoms

    return

end
