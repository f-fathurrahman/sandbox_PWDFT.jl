function calc_QVeff!( Ham )

    Nspin = Ham.electrons.Nspin
    Ng = Ham.pw.gvec.Ng
    lmaxq = Ham.pspotNL.lmaxq
    G2 = Ham.pw.gvec.G2
    G = Ham.pw.gvec.G
    Npoints = prod(Ham.pw.Ns)
    idx_g2r = Ham.pw.idx_g2r
    Veff = Ham.potentials.Total

    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    # Ylm_real_qe accept l value starting from 0
    Ylm_real_qe!(_lmax, G, ylmk0)

    VeffG = zeros(ComplexF64, Ng, Nspin)

    # Fourier transform of the total effective potential
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        for ip in 1:Npoints
            ctmp[ip] = Veff[ip,ispin] # Veff already contains Ps_loc
        end
        R_to_G!(Ham.pw, ctmp)
        end ig in 1:Ng
            ip = idx_g2r[ig]
            VeffG[ig,ispin] = ctmp[ip]
        end
    end


    Nspecies = Ham.atoms.Nspecies
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    nh = Ham.pspotNL.nh
    Deeq = Ham.pspotNL.Deeq
    CellVolume = Ham.pw.CellVolume

    for isp in 1:Nspecies

        if Ham.pspots[isp].is_ultrasoft
            # nij = max number of (ih,jh) pairs per atom type nt
            nij = nh[isp] * (nh[isp] + 1)/2

            Qgm = zeros(ComplexF64, Ng, nij)
            #
            # Compute and store Q(G) for this atomic species 
            # (without structure factor)
            ijh = 0
            for ih in 1:nh[isp], jh in ih:nh[isp]
                ijh = ijh + 1
                #qvan2!( ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
                @views qvan2!( Ham.pspotNL, ih, jh, isp, G2, ylmk0, Qgm[:,ijh] )
            end
            #
            # count max number of atoms of type nt
            #
            nab = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                nab = nab + 1
            end
            #
            aux = zeros(ComplexF64, Ng, nab)
            #
            # Compute and store V(G) times the structure factor e^(-iG*tau)
            for ispin in 1:Nspin
                nb = 0
                for ia in 1:Natoms
                    if atm2species[ia] != isp
                        continue
                    end
                    nb = nb + 1
                    for ig in 1:Ng
                        GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
                        Sf = cos(GX) - im*sin(GX)
                        aux[ig,nb] = VeffG[ig,ispin] * Sf
                    end
                end
            #
            # here we compute the integral Q*V for all atoms of this kind
            deeaux = Qgm*aux
            #CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, aux, &
            #         2*ngm, 0.0_dp, deeaux, nij )
        
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                nb = nb + 1
                ijh = 0
                for ih in 1:nh[isp], jh = ih:nh[isp]
                    ijh = ijh + 1
                    Deeq[ih,jh,ia,ispin] = CellVolume * deeaux[ijh,nb]
                    if jh > ih
                        Deeq[jh,ih,ia,ispin] = Deeq[ih,jh,ia,ispin]
                    end
                end
            end # Natoms
        
        end # is_ultrasoft
    
    end # Nspecies

end