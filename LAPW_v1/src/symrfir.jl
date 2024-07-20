function symrfir!(pw, sym_vars, rfir)

    Npoints = prod(pw.Ns)
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    idx_g2miller = pw.gvec.idx_g2miller

    nsymcrys = sym_vars.nsymcrys
    vtcsymc = sym_vars.vtcsymc
    tv0symc = sym_vars.tv0symc
    lsplsymc = sym_vars.lsplsymc
    symlat = sym_vars.symlat
    isymlat = sym_vars.isymlat

    # Fourier transform function to G-space
    zfft1 = zeros(ComplexF64, Npoints)
    zfft1[:] = rfir[:]
    R_to_G!(pw, zfft1)

    zfft2 = zeros(ComplexF64, Npoints)
    # loop over crystal symmetries
    for isym in 1:nsymcrys
        # zero translation vector flag
        tv0 = tv0symc[isym]
        # translation vector in Cartesian coordinates
        if !tv0 
            v1 = vtcsymc[1,isym]
            v2 = vtcsymc[2,isym]
            v3 = vtcsymc[3,isym]
        end
        # index to lattice symmetry of spatial rotation
        lspl = lsplsymc[isym]
        # inverse rotation required for rotation of G-vectors
        ilspl = isymlat[lspl]
        sym = symlat[ilspl]
        # Begin loop over all G-vectors
        for ig in 1:Ng
            ip = idx_g2r[ig]
            # multiply the transpose of the inverse symmetry matrix with the G-vector
            if lspl == 1
                jp = ip
                # No need to rotate ?
            else
                i1 = idx_g2miller[ig][1]
                i2 = idx_g2miller[ig][2]
                i3 = idx_g2miller[ig][3]
                # multiplication by inverse or transpose of rotation matrix
                j1 = sym[1,1]*i1 + sym[2,1]*i2 + sym[3,1]*i3
                j2 = sym[1,2]*i1 + sym[2,2]*i2 + sym[3,2]*i3
                j3 = sym[1,3]*i1 + sym[2,3]*i2 + sym[3,3]*i3
                # Find idx of the rotated G vector
                # XXX: We also can use G2_shells here
                jg = findfirst( isequal((j1,j2,j3)), idx_g2miller )
                jp = idx_g2r[jg]
            end
            if tv0 
                # zero translation vector
                zfft2[jp] += zfft1[ip]
            else
                # complex phase factor for translation
                t1 = -( G[1,ig]*v1 + G[2,ig]*v2 + G[3,ig]*v3 )
                # XXX: Compare with QE
                z1 = cos(t1) + im*sin(t1)
                zfft2[jp] += z1*zfft1[ip]
            end 
        end
    end
    # Fourier transform to real-space
    G_to_R!(pw, zfft2)
    # Normalize by number of crystal symmetries
    t1 = 1.0/nsymcrys
    rfir[:] = t1*real(zfft2)

    return
end # function 

