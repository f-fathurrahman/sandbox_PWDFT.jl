# XXX implement findsymlat here directly
function findsymlat!(atoms, sym_vars)
    nsymlat, symlat, symlatc, symlatd, isymlat = findsymlat(atoms)
    sym_vars.nsymlat = nsymlat
    sym_vars.symlat = symlat
    sym_vars.symlatc = symlatc
    sym_vars.symlatd = symlatd
    sym_vars.isymlat = isymlat
    return
end



#  Finds the point group symmetries which leave the Bravais lattice invariant.
#  Let $A$ be the matrix consisting of the lattice vectors in columns, THEN 
#  $$ g=A^{\rm T}A $$
#  is the metric tensor. Any $3\times 3$ matrix $S$ with elements $-1$, 0 or 1
#  is a point group symmetry of the lattice if $\det(S)$ is $-1$ or 1, and
#  $$ S^{\rm T}gS=g. $$
#  The first matrix in the set RETURN ed is the identity.
function findsymlat(atoms; epslat=1e-6)
    # XXX: only pass LatVecs instead of atoms

    LatVecs = atoms.LatVecs
    # Metric tensor
    g = LatVecs' * LatVecs
    
    #=
    INTEGER :: md, sym(3,3), i, j
    INTEGER :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    REAL(8) :: s(3,3), g(3,3), sgs(3,3)
    REAL(8) :: c(3,3), v(3), t1
    ! external functions
    INTEGER :: i3mdet
    external i3mdet
    =#

    MAX_SYM_LATT = 48

    S = zeros(Int64, 3, 3)
    SgS = zeros(Float64, 3, 3)
    symlat = Vector{Matrix{Int64}}(undef,MAX_SYM_LATT)
    for i in 1:48
        symlat[i] = zeros(Int64, 3, 3)
    end
    symlatd = zeros(Int64, MAX_SYM_LATT)
    isymlat = zeros(Int64, MAX_SYM_LATT)
    #
    # loop over all possible symmetry matrices
    nsymlat = 0
    for i11 in -1:1, i12 in -1:1, i13 in -1:1
        for i21 in -1:1, i22 in -1:1, i23 in -1:1
            for i31 in -1:1, i32 in -1:1, i33 in -1:1
                #
                S[1,1] = i11; S[1,2] = i12; S[1,3] = i13
                S[2,1] = i21; S[2,2] = i22; S[2,3] = i23
                S[3,1] = i31; S[3,2] = i32; S[3,3] = i33
                # determinant of matrix
                md = Int64(det(S)) # convert to Int64
                #
                # matrix should be orthogonal
                if abs(md) != 1
                    continue
                end
                #
                # check invariance of metric tensor
                @views SgS[:,:] = (S' * g) * S
                if any( abs.(SgS .- g) .> epslat )
                    continue
                end 
                #
                # additional checks for spin-spiral, electric field, and A-field is skipped
                # ....
                #
                nsymlat += 1
                if nsymlat > 48
                    error("Too many nsymlat: nsymlat = $(nsymlat)")
                end 
                # 
                symlat[nsymlat][:,:] = S[:,:]
                symlatd[nsymlat] = md
            end
        end
    end
    
    @info "After loop_matrix33 nsymlat = $nsymlat"
    for i in 1:nsymlat
        @info "$i symlat = $(symlat[i])"
    end

    if nsymlat == 0 
        error("Lattice symmetry not found")
    end
    
    # make the first symmetry the identity
    for i in 1:nsymlat
        if symlat[i] == I(3)
            S[:,:] = symlat[1]
            symlat[1][:,:] = symlat[i][:,:]
            symlat[i][:,:] = S[:,:]
            md = symlatd[1]
            symlatd[1] = symlatd[i]
            symlatd[i] = md
            break
        end
    end
    #
    # index to the inverse of each operation
    # XXX: Probably it is simpler to just compute the inverse
    #      This might be useful to check if symlat elements form a group
    for i in 1:nsymlat
        @views S[:,:] = Int64.(inv(symlat[i]))
        for j in 1:nsymlat
            if S == symlat[j]
                isymlat[i] = j
                @goto INVERSE_FOUND_LABEL
            end
        end 
        error("Inverse operation is not found")
        @label INVERSE_FOUND_LABEL
    end
    #
    # determine the lattice symmetries in Cartesian coordinates
    symlatc = Vector{Matrix{Float64}}(undef,nsymlat)
    for i in 1:nsymlat
        symlatc[i] = zeros(Float64, 3, 3)
    end
    ainv = inv(LatVecs)
    for i in 1:nsymlat
        symlatc[i] = LatVecs * ( symlat[i] * ainv )
    end
    return nsymlat, symlat[1:nsymlat], symlatc, symlatd[1:nsymlat], isymlat[1:nsymlat]
end

