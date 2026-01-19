function symmetrize_rhoe_v2!(
    pw::PWGrid,
    sym_info::SymmetryInfo,
    rhoe_symmetrizer::RhoeSymmetrizer,
    Rhoe::Array{Float64,2}
)

    Ngs = rhoe_symmetrizer.Ngs
    shell_G_sym = rhoe_symmetrizer.shell_G_sym
    fts = rhoe_symmetrizer.fts

    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    Nsyms = sym_info.Nsyms
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    non_symmorphic = sym_info.non_symmorphic

    Npoints = prod(pw.Ns)
    Nspin = size(Rhoe)[2]

    RhoeG = zeros(ComplexF64,Npoints,Nspin)
    for ispin = 1:Nspin
        @views RhoeG[:,ispin] = R_to_G(pw, Rhoe[:,ispin])
    end

    sg = zeros(Float64,3)
    irot = zeros(Int64,Nsyms)

    g0 = zeros(3,Nsyms)
    is_done_shell = zeros(Bool,Nsyms)

    trmat = LatVecs'/(2*pi)

    rhosum = zeros(ComplexF64,Nspin)
    
    non_symmorphic = zeros(Bool, 48)
    #
    # convert fractional translations to cartesian, in a0 units
    for isym in 1:Nsyms
        non_symmorphic[isym] = ( ft[1,isym] != 0.0 || ft[2,isym] != 0.0 || ft[3,isym] != 0.0 )
        if non_symmorphic[isym]
            ft_[:,isym] = at[:,1]*ft[1,isym] + at[:,2]*ft(2,ns) + at[:,3]*ft[3,isym]
        end
    end
    
    if Nspin_dens == 4
        Nspin_lsda = 1
    elseif Nspin_dens in [1,2]
        Nspin_lsda = Nspin_dens
    else
        @error("Not valid number of Nspin_dens=$(Nspin_dens)")
    end
    #
    # scan shells of G-vectors
    #
    for igl in 1:Ngs
        
        # symmetrize: \rho_sym(G) = \sum_S rho(SG) for all G-vectors in the star
        #
        ng = length(shell_G_sym[igl])
        @assert ng >= 1
        for ig in 1:ng
            idx_shell = shell_G_sym[igl][ig]
            g0[1,ig] = G[1,idx_shell]
            g0[2,ig] = G[2,idx_shell]
            g0[3,ig] = G[3,idx_shell]
            is_done_shell[ig] = false
        end
        g0[:,:] = trmat*g0[:,:]
        #
        #  rotate G-vectors
        #
        for ig in 1:ng
            #
            if is_done_shell[ig]
                continue
            end
            
            fill!(rhosum, 0.0 + im*0.0)
            fill!(magsum, 0.0 + im*0.0)

            for isym in 1:Nsyms
                for ii in 1:3
                    sg[ii] = inv_s[ii,1,isym]*g0[1,ig] + inv_s[ii,2,isym]*g0[2,ig] + inv_s[ii,3,isym]*g0[3,ig]
                end
                #
                irot[isym] = 0
                for isg in 1:ng
                    if (abs(sg[1] - g0[1,isg]) < 1e-5) &&
                       (abs(sg[2] - g0[2,isg]) < 1e-5) &&
                       (abs(sg[3] - g0[3,isg]) < 1e-5)
                        irot[isym] = isg
                        break
                    end
                end

                if (irot[isym] < 1) || (irot[isym] > ng)
                    error("Error in determining irot")
                end

                # isg is the index of rotated G-vector
                isg = shell_G_sym[igl][irot[isym]]
                ip = idx_g2r[isg]
                #
                # non-spin-polarized case: component 1 is the charge
                # LSDA case: components 1,2 are spin-up and spin-down charge
                # non colinear case: component  1 is the charge density,
                #                    components 2,3,4 are the magnetization
                # non colinear case: components 2,3,4 are the magnetization
                #
                if Nspin_dense == 4
                    # bring magnetization to crystal axis
                    mag[:] = RhoeG[ip,2] * RecVecs[1,:] +
                             RhoeG[ip,3] * RecVecs[2,:] +
                             RhoeG[ip,4] * RecVecs[3,:]
                    # rotate and add magnetization
                    magrot[:] = s[1,:,invs[isym]] * mag[1] +
                                s[2,:,invs[isym]] * mag[2] +
                                s[3,:,invs[isym]] * mag[3]
                    if sname[invs[isym]][1:3] == "inv"
                        magrot[:] = -magrot[:]
                    end
                    if t_rev[invs[isym]]
                        magrot[:] = -magrot[:]
                    end
                end

                if non_symmorphic[isym]
                    arg = G[1,isg]*fts[1,isym] + G[2,isg]*fts[2,isym] + G[3,isg]*fts[3,isym]
                    fact = cos(arg) - im*sin(arg)
                    for ispin in 1:Nspin_lsda
                        rhosum[ispin] += RhoeG[ip,ispin]*fact
                    end
                    if Nspin_dens == 4
                        magsum[:] .+= magrot[:] * fact
                    end
                else
                    for ispin in 1:Nspin_lsda
                        rhosum[ispin] += RhoeG[ip,ispin]
                    end
                    if Nspin_dens == 4
                        magsum[:] .+= magrot[:]
                    end
                end
                
            end # isym
            #
            for ispin in 1:Nspin_lsda
                rhosum[ispin] = rhosum[ispin]/Nsyms
            end

            for isym in 1:Nsyms
                isg = shell_G_sym[igl][irot[isym]]
                ip = idx_g2r[isg]
                if non_symmorphic[isym]
                    arg = G[1,isg]*fts[1,isym] + G[2,isg]*fts[2,isym] + G[3,isg]*fts[3,isym]
                    fact = cos(arg) + im*sin(arg)
                    for ispin in 1:Nspin_lsda
                        RhoeG[ip,ispin] = rhosum[ispin]*fact
                    end
                else
                    for ispin in 1:Nspin
                        RhoeG[ip,ispin] = rhosum[ispin]
                    end
                end
                is_done_shell[irot[isym]] = true
            end

        end # ng
    end # ngl

    for ispin in 1:Nspin_dens
        @views Rhoe[:,ispin] = real(G_to_R(pw, RhoeG[:,ispin]))
    end

    return
end
