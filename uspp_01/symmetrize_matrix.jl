function symmetrize_matrix!( LatVecs_, sym_info::SymmetryInfo, M )
    # only for 3x3 matrix
    # Symmetrize a function \(f(i,j)\) (e.g. stress, dielectric tensor in
    # cartesian axis), where \(i,j\) are the cartesian components.
    
    Nsyms = sym_info.Nsyms
    if Nsyms == 1
        return
    end

    alat = norm(LatVecs_[:,1])
    LatVecs = LatVecs_[:,:]/alat
    RecVecs = inv(Matrix(LatVecs'))

    # bring matrix to crystal axis
    # CALL cart_to_crys( matr )
    M_cryst = zeros(Float64, 3, 3)
    for i in 1:3
        M_cryst[:,i] = M[1,i]*LatVecs[1,:] + M[2,i]*LatVecs[2,:] + M[3,i]*LatVecs[3,:]
    end
    # XXX: Check this!

    s = convert(Array{Float64,3}, sym_info.s)
    # symmetrize in crystal axis
    work = zeros(Float64, 3, 3)
    for isym in 1:Nsyms
        for i in 1:3, j in 1:3, k in 1:3, for l in 1:3
            work[i,j] += s[i,k,isym] * s[j,l,isym] * M_cryst[k,l]
        end
    end
    M_cryst[:,:] = work[:,:] / Nsyms
    
    # bring matrix back to cartesian axis
    for i in 1:3
        M[:,i] = M_cryst[1,i]*RecVecs[:,1] + M_cryst[2,i]*RecVecs[:,2] + M_cryst[3,i]*RecVecs[:,3]
    end

    return
end