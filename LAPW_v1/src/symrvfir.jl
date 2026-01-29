function symrvfir!(pw, sym_vars, rvfir; tspin = true, tnc = false)

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
    lspnsymc = sym_vars.lspnsymc
    symlatd = sym_vars.symlatd
    symlatc = sym_vars.symlatc

    # dimension of the vector field
    if tnc
      nd = 3
    else
      nd = 1
    end
  
    zfft1 = zeros(ComplexF64, Npoints, nd)
    zfft2 = zeros(ComplexF64, Npoints, nd)

    sc = zeros(Float64, 3, 3)
    zv1 = zeros(ComplexF64, 3)
    zv2 = zeros(ComplexF64, 3)
    sym = zeros(Int64, 3, 3)

    # Fourier transform vector function to G-space
    for i in 1:nd
        zfft1[:,i] = rvfir[:,i]
        @views R_to_G!(pw, zfft1[:,i])
    end

    for isym in 1:nsymcrys
        # zero translation vector flag
        tv0 = tv0symc[isym]
        # translation vector in Cartesian coordinates
        if !tv0
            v1 = vtcsymc[1,isym]
            v2 = vtcsymc[2,isym]
            v3 = vtcsymc[3,isym]
        end
        # index to spatial rotation lattice symmetry
        lspl = lsplsymc[isym]
        # inverse rotation required for rotation of G-vectors
        ilspl = isymlat[lspl]
        sym[:,:] = symlat[ilspl]
        if tspin
            # global spin proper rotation in Cartesian coordinates
            lspn = lspnsymc[isym]
            sc[:,:] = symlatd[lspn]*symlatc[lspn]
        else
            # set spin rotation equal to spatial rotation
            lspn = lspl
            sc[:,:] = symlatc[lspl]
        end
        for ig in 1:Ng
            ip = idx_g2r[ig]
            # multiply the transpose of the inverse symmetry matrix with the G-vector
            if lspl == 1
                jp = ip
            else
                i1, i2, i3 = idx_g2miller[ig]
                j1 = sym[1,1]*i1 + sym[2,1]*i2 + sym[3,1]*i3
                j2 = sym[1,2]*i1 + sym[2,2]*i2 + sym[3,2]*i3
                j3 = sym[1,3]*i1 + sym[2,3]*i2 + sym[3,3]*i3
                jg = findfirst( isequal((j1,j2,j3)), idx_g2miller )
                jp = idx_g2r[jg]
            end
            # complex phase factor for translation
            if tv0
                z1 = 1.0
            else
                t1 = -(G[1,ig]*v1 + G[2,ig]*v2 + G[3,ig]*v3)
                z1 = cos(t1) + im*sin(t1)
            end
            # translation, spatial rotation and global spin rotation
            if lspn == 1
                # global spin symmetry is the identity
                zfft2[jp,:] += z1*zfft1[ip,:]
            else
                if tnc
                    # non-collinear case
                    zv1[:] = zfft1[ip,:]
                    zv2[1] = sc[1,1]*zv1[1] + sc[1,2]*zv1[2] + sc[1,3]*zv1[3]
                    zv2[2] = sc[2,1]*zv1[1] + sc[2,2]*zv1[2] + sc[2,3]*zv1[3]
                    zv2[3] = sc[3,1]*zv1[1] + sc[3,2]*zv1[2] + sc[3,3]*zv1[3]
                    zfft2[jp,:] = zfft2[jp,:] + z1*zv2[:]
                else
                    # collinear case
                    zfft2[jp,1] += sc[3,3]*z1*zfft1[ip,1]
                end
            end
        end
    end
    # Fourier transform to real-space and normalize
    t1 = 1.0/nsymcrys
    for i in 1:nd
        @views G_to_R!(pw, zfft2[:,i])
        rvfir[:,i] = t1*real.(zfft2[:,i])
    end
    return
end
