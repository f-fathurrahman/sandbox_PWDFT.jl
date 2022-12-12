# From d_matrix.f90 in QE

# sr: rotation matrices in Cartesian
function init_d_matrices( sr::Array{Float64,3} )
    # !! Calculates the d-matrices.
    # !
    # USE kinds,            ONLY: DP
    # USE symm_base,        ONLY: nsym, sr
    # USE random_numbers,   ONLY: randy
    # USE matrix_inversion

    Nsyms = size(sr, 3)

    # Hardcoded parameters 
    # MAXL = max value of l allowed
    # MAXL = number of m components for l=MAXL
    # MAXLM = number of l,m spherical harmonics for l <= MAXL
    MAXL = 3
    MAXM = 2*MAXL + 1
    MAXLM = (MAXL + 1)^2

    SMALL = 1.0e-9

    # These arrays will be returned
    dy1 = zeros(Float64, 3, 3, 48)
    dy2 = zeros(Float64, 5, 5, 48)
    dy3 = zeros(Float64, 7, 7, 48)

    #
    # randomly distributed points on a sphere
    #
    rl = zeros(Float64, 3, MAXM)
    for m in 1:MAXM
        rl[1,m] = rand() - 0.5
        rl[2,m] = rand() - 0.5
        rl[3,m] = rand() - 0.5
    end
  
    # CALL ylmr2( maxlm, 2*maxl+1, rl, rrl, ylm )
    ylm = zeros(Float64, MAXM, MAXLM)
    Ylm_real_qe!(MAXL, rl, ylm) # Ylm_real_qe accept l value starting from 0


    # invert Yl for each block of definite l (note the transpose operation)
    
    # l = 1 block
    yl1 = zeros(Float64, 3, 3)
    for m in 1:3, n in 1:3
        yl1[m,n] = ylm[n,1+m]
    end
    yl1_inv = inv(yl1)  
    
    #  l = 2 block
    yl2 = zeros(Float64, 5, 5)
    for m in 1:5, n in 1:5
        yl2[m,n] = ylm[n,4+m]
    end
    yl2_inv = inv(yl2)


    #  l = 3 block
    yl3 = zeros(Float64, 7, 7)
    for m in 1:7, n in 1:7
        yl3[m,n] = ylm[n,9+m]
    end
    yl3_inv = inv(yl3)



    # now for each symmetry operation of the point-group

    srl = zeros(Float64, 3, MAXM)
    delt = zeros(Float64, 7, 7)
    ylms = zeros(Float64, MAXM, MAXLM)

    for isym in 1:Nsyms
        #
        # srl[:,m] = rotated rl[:,m] vectors
        #
        @views srl[:,:] = sr[:,:,isym] * rl
        # srl = MATMUL( sr(:,:,isym), rl )
        
        # CALL ylmr2( maxlm, maxm, srl, rrl, ylms )
        Ylm_real_qe!(MAXL, srl, ylms)

        #  find  D_S = Yl_S * Yl_inv (again, beware the transpose)
        
        # l = 1
        for m in 1:3, n in 1:3
            yl1[m,n] = ylms[n,1+m]
        end
        @views dy1[:,:,isym] = yl1[:,:] * yl1_inv[:,:]

        # l = 2 block
        for m in 1:5, n in 1:5
            yl2[m,n] = ylms[n,4+m]
        end
        @views dy2[:,:,isym] = yl2[:,:] * yl2_inv[:,:]


        # l = 3 block
        for m in 1:7, n in 1:7
            yl3[m,n] = ylms[n,9+m]
        end
        @views dy3[:,:,isym] = yl3[:,:] * yl3_inv[:,:]
    end


    # check that D_S matrices are orthogonal as they should if Ylm are correctly defined.
    fill!( delt, 0.0 )
    for m in 1:7
        delt[m,m] = 1.0
    end
    
    for isym in 1:Nsyms
        _check_is_d_orthogonal(isym, dy1, 1, 3, delt, SMALL)
        _check_is_d_orthogonal(isym, dy2, 2, 5, delt, SMALL)
        _check_is_d_orthogonal(isym, dy3, 3, 7, delt, SMALL)
    end
    return dy1, dy2, dy3
end


function _check_is_d_orthogonal(isym, D, l, Nsize, delt, SMALL)
    # l = 3 block
    cc = 0.0
    for m in 1:Nsize, n in 1:Nsize
        dd = dot(D[:,m,isym], D[:,n,isym])
        cc += ( dd - delt[m,n] )^2
    end
    if cc > SMALL
        @printf("D_S (l=%d) for this symmetry operation (isym=%d) is not orthogonal\n", l, isym)
        error("Error in init_d_matrices")
    end
    return
end
