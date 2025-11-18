function symrvfmt!(tspin,tnc,nr,nri,np,ld,rvfmt)
    
    tnc = false
    # dimension of the vector field
    if tnc
        nd = 3
    else
        nd = 1
    end

    rvfmt1 = zeros(Float64, npmtmax, Natoms, nd)
    rvfmt2 = zeros(Float64, npmtmax, nd)
    sc = zeros(Float64, 3, 3)
    v1 = zeros(Float64, 3)
    v2 = zeros(Float64, 3)
    #
    t0 = 1.0/nsymcrys
    for isp in 1:nspecies
        # make copy of vector field for all atoms of current species
        for i in 1:nd, ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            rvfmt1[:,ia,i] = rvfmt[ia][:,i]
        end
        done[:] = false
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            #
            if done[ia]
                continue
            end
            rvfmt[ia][1:np[isp],1:nd] .= 0.0
            # begin loop over crystal symmetries
            for isym in 1:nsymcrys
                # equivalent atom
                ja = ieqatom[ia][isym]
                # parallel transport of vector field
                lspl = lsplsymc[isym]
                for i in 1:nd
                    @views rotrfmt!(mt_vars, isp, symlatc[lspl], rvfmt1[:,ja,i], rvfmt2[:,i])
                end
                #
                if tspin
                    # global spin proper rotation matrix in Cartesian coordinates
                    lspn = lspnsymc[isym]
                    sc[:,:] = symlatd[lspn]*symlatc[lspn]
                else
                    # set spin rotation equal to spatial rotation
                    lspn = lspl
                    sc[:,:] = symlatc[lspls]
                end
                # global spin rotation of vector field
                if tnc
                    # non-collinear case
                    for i in 1:np[isp]
                        v1[:] = rvfmt2[i,:]
                        v2[1] = sc[1,1]*v1[1] + sc[1,2]*v1[2] + sc[1,3]*v1[3]
                        v2[2] = sc[2,1]*v1[1] + sc[2,2]*v1[2] + sc[2,3]*v1[3]
                        v2[3] = sc[3,1]*v1[1] + sc[3,2]*v1[2] + sc[3,3]*v1[3]
                        @. rvfmt[ia][i,1:3] += v2[1:3]
                    end
                else
                    # collinear case
                    #call daxpy(np(is),sc(3,3),rvfmt2,1,rvfmt(:,ias,1),1)
                    rvfmt[ia][:,ias] += sc[3,3]*rvfmt2[1:np[isp]]
                end
            end # end loop over crystal symmetries
            # normalize
            for i in 1:nd
                #call dscal(np(is),t0,rvfmt(:,ias,i),1)
                rvfmt[ia][:,i] *= t0
            end
            # mark atom as done
            done[ia] = true
            # rotate into equivalent atoms
            for isym in 1:nsymcrys
                ja = ieqatom[ia][isym]
                if done[ja]
                    continue
                end
                # parallel transport of vector field (using operation inverse)
                lspl = isymlat(lsplsymc(isym))
                for i in 1:nd
                    @views rotrfmt!(mt_vars, isp, symlatc[lspl], rvfmt[ia][:,i], rvfmt[ja][:,i])
                end
                if tspin
                    # inverse of global proper rotation matrix in Cartesian coordinates
                    lspn = isymlat[lspnsymc[isym]]
                    sc[:,:] = symlatd[lspn]*symlatc[lspn]
                else
                    # set spin rotation equal to spatial rotation
                    lspn = lspl
                    sc[:,:] = symlatc[lspl]
                end
                # global spin rotation of vector field
                if tnc
                    # non-collinear case
                    for i in 1:np[is]
                        v1[1:3] = rvfmt[ja][i,1:3]
                        v2[1] = sc[1,1]*v1[1] + sc[1,2]*v1[2] + sc[1,3]*v1[3]
                        v2[2] = sc[2,1]*v1[1] + sc[2,2]*v1[2] + sc[2,3]*v1[3]
                        v2[3] = sc[3,1]*v1[1] + sc[3,2]*v1[2] + sc[3,3]*v1[3]
                        rvfmt[ja][i,1:3] = v2[1:3]
                    end
                else
                    # collinear case
                    #call dscal(np(is),sc(3,3),rvfmt(:,jas,1),1)
                    rvfmt[ja][:,1] *= sc[3,3]
                end
                # mark atom as done
                done[ja] = true
            end
        end  # end loop over atoms and species
    end
    return
end
