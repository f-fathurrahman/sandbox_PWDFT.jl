#=
Finds the symmetries which rotate one set of atomic positions into another.
Both sets of positions differ only by a translation vector and have the same
muffin-tin magnetic fields (stored in the global array {\tt bfcmt}). Any
symmetry element consists of a spatial rotation of the atomic position
vectors followed by a global magnetic rotation: $\{\alpha_S|\alpha_R\}$. In
the case of spin-orbit coupling $\alpha_S=\alpha_R$. The symmetries are
returned as indices of elements in the Bravais lattice point group. An
index to equivalent atoms is stored in the array {\tt iea}.
=#
function findsym!(
    sym_vars::SymmetryVars,
    atoms::Atoms,
    apl1, apl2,
    lspl, lspn, iea;
    epslat = 1e-6,
    spinpol = false,
    spinorb = false,
    bfcmt0::Union{Matrix{Float64}, Nothing} = nothing,
    bfieldc0::Union{Vector{Float64}, Nothing} = nothing
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nsymlat = sym_vars.nsymlat
    symlat = sym_vars.symlat
#=
    println()
    println("--- ENTER findsym ---")
    println()

    println()
    println("First positions")
    for ia in 1:Natoms
        @printf("%18.10f %18.10f %18.10f\n", apl1[1,ia], apl1[2,ia], apl1[3,ia])
    end

    println()
    println("Second positions")
    for ia in 1:Natoms
        @printf("%18.10f %18.10f %18.10f\n", apl2[1,ia], apl2[2,ia], apl2[3,ia])
    end

    println()
=#

    jea = zeros(Int64, Natoms)
    sl = zeros(Float64, 3, 3)
    apl3 = zeros(Float64, 3, Natoms)
    v = zeros(Float64, 3)

    # These variables are not used if spinpol=true
    sc = zeros(Float64, 3, 3)
    if isnothing(bfcmt0)
        bfcmt0 = zeros(Float64, 3, Natoms)
    end
    if isnothing(bfieldc0)
        bfieldc0 = zeros(Float64, 3)
    end


    nsym = 0
    sl = zeros(Float64, 3, 3)
    # loop over lattice symmetries (spatial rotations)
    for isym in 1:nsymlat
        #
        #println("\nisym = ", isym)
        #
        # make real copy of lattice rotation symmetry
        @views sl[:,:] = Float64.(symlat[isym][:,:])
        # loop over species
        for ia in 1:Natoms
            isp = atm2species[ia]
            # map apl1 coordinates to [0,1) and store in apl3
            @views apl3[:,ia] = apl1[:,ia]
            @views r3frac!(apl3[:,ia], epslat=epslat)
            #
            for ja in 1:Natoms
                if atm2species[ja] != isp
                    #@info "Skipping this atom index"
                    continue
                end
                # apply lattice symmetry to atomic positions
                @views v[:] = sl[:,1]*apl2[1,ja] + sl[:,2]*apl2[2,ja] + sl[:,3]*apl2[3,ja]
                # map coordinates to [0,1)
                r3frac!(v, epslat=epslat)
                # check if atomic positions are invariant
                t1 = abs(apl3[1,ia]-v[1]) + abs(apl3[2,ia]-v[2]) + abs(apl3[3,ia]-v[3])
                #@info "t1 = $t1"
                #println("Checking isp=$isp ia=$ia ja=$ja")
                if t1 < epslat 
                    #println(" *** Equivalent atoms: ia=$ia ja=$ja isym=$isym")
                    # equivalent atom index
                    jea[ia] = ja
                    @goto LABEL10 # continue ?
                end
            end
            # not invariant so try new spatial rotation
            #println(" - Not invariant, trying new spatial rotation")
            @goto LABEL40 # jump to the end of loop over symlat
            @label LABEL10 #10 CONTINUE
        end
    
        # all atomic positions invariant at this point
        jsym = 1

        # 
        # .... spin polarized stuffs is removed
        # spinpol
        if spinpol
            symlatd = sym_vars.symlatd
            symlatc = sym_vars.symlatc
            # check invariance of magnetic fields under global spin rotation
            if spinorb 
                # with spin-orbit coupling spin rotation equals spatial rotation
                jsym0 = isym
                jsym1 = isym
            else
                # without spin-orbit coupling spin rotation independent of spatial rotation
                jsym0 = 1
                jsym1 = nsymlat
            end
            #
            for jsym in jsym0:jsym1
                # determinant of the symmetry matrix
                md = symlatd[jsym]
                sc[:,:] = md*symlatc[jsym]
                # rotate global field and check invariance using proper part of symmetry matrix
                v[:] = sc[:,1]*bfieldc0[1] + sc[:,2]*bfieldc0[2] + sc[:,3]*bfieldc0[3]
                t1 = abs(bfieldc0[1]-v[1]) + abs(bfieldc0[2] - v[2]) + abs(bfieldc0[3] - v[3])
                # if not invariant try a different global spin rotation
                if t1 > epslat
                    @goto LABEL20
                end
                # rotate muffin-tin magnetic fields and check invariance
                for ia in 1:Natoms
                    isp = atm2species[ia]
                    # equivalent atom
                    ja = jea[ia]
                    v[:] = sc[:,1]*bfcmt0[1,ja] + sc[:,2]*bfcmt0[2,ja] + sc[:,3]*bfcmt0[3,ja]
                    t1 = abs(bfcmt0[1,ia] - v[1]) + abs(bfcmt0[2,ia] - v[2]) + abs(bfcmt0[3,ia] - v[3])
                    # if not invariant try a different global spin rotation
                    if t1 > epslat
                        @goto LABEL20
                    end 
                end
                # all fields invariant
                @goto LABEL30
                @label LABEL20 # continue
            end # end loop over global spin rotations 
            # magnetic fields not invariant so try different spatial rotation
            @goto LABEL40  # nsym is not incremented
        end

        @label LABEL30
        # everything invariant so add symmetry to set
        nsym += 1
        lspl[nsym] = isym
        lspn[nsym] = jsym
        for ia in 1:Natoms
            iea[ia,nsym] = jea[ia]
        end
        #println("nsym = ", nsym)

        #@info "Trying new spatial rotation"
        @label LABEL40 # continue
    
    end # ! end loop over spatial rotations 

    #println("nsym = ", nsym)

    #println()
    #println("--- EXIT findsym ---")
    #println()

    return nsym

    return
end

