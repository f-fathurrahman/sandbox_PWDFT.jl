include("calc_Deff.jl")


function my_calc_forces_Ps_nloc!(
    atoms::Atoms,    
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    electrons::Electrons,
    pspotNL::PsPotNL_UPF,
    potentials::Potentials,
    psiks::BlochWavefunc,
    F_Ps_nloc::Matrix{Float64}
)

    betaNL = pspotNL.betaNL
    idx_gw2g = pw.gvecw.idx_gw2g
    G = pw.gvec.G
    NbetaNL = pspotNL.NbetaNL
    Nstates = size(psiks[1],2) # or electrons.Nstates
    Nspin = electrons.Nspin
    Nkpt = pw.gvecw.kpoints.Nkpt

    Ngw = pw.gvecw.Ngw
    Ngwx = maximum(Ngw)
    # To save memory usage
    dbetaNL = zeros(ComplexF64,Ngwx,NbetaNL)

    betaNL_psi = zeros(ComplexF64,NbetaNL,Nstates)
    dbetaNL_psi = zeros(ComplexF64,NbetaNL,Nstates)

    fill!(F_Ps_nloc, 0.0)

    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        Ngw_ik = Ngw[ik]
        #
        betaNL_psi[:,:] .= betaNL[ik]' * psi
        # Calculate derivative of betaNL in G-space
        for ipol in 1:3
            for ibeta in 1:NbetaNL, igw in 1:Ngw_ik
                ig = idx_gw2g[ik][igw]
                dbetaNL[igw,ibeta] = -im * betaNL[ik][igw,ibeta] * G[ipol,ig]
            end
            # betaNL psi 
            @views dbetaNL_psi[:,:] .= dbetaNL[1:Ngw_ik,:]' * psi
            #
            # this will call sum over bands
            _force_Ps_nloc_k!(ipol, ik, ispin, atoms, pw, pspots,
                electrons, pspotNL, betaNL_psi, dbetaNL_psi, F_Ps_nloc)
        end
    end

    # This does not depend on k-points
    _add_F_uspp!(atoms, pw, pspots, pspotNL, potentials, F_Ps_nloc)


    return

end


function _force_Ps_nloc_k!(ipol, ik, ispin, 
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    electrons::Electrons,
    pspotNL::PsPotNL_UPF,
    betaNL_psi, dbetaNL_psi,
    F_Ps_nloc
)

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies

    ebands = electrons.ebands
    Nstates = electrons.Nstates
    Focc = electrons.Focc

    nh = pspotNL.nh
    nhm = pspotNL.nhm
    Deeq = pspotNL.Deeq
    indv_ijkb0 = pspotNL.indv_ijkb0

    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk

    Deff = zeros(Float64, nhm, nhm, Natoms)

    ikspin = ik + (ispin - 1)*Nkpt
    #
    for ist in 1:Nstates
        #
        _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ik], Deff )
        #
        fac = wk[ik] * Focc[ist,ikspin]
        #
        for isp in 1:Nspecies, ia in 1:Natoms
            
            if atm2species[ia] != isp
                continue
            end

            psp = pspots[isp]            
            ijkb0 = indv_ijkb0[ia]

            for ih in 1:nh[isp]
                ikb = ijkb0 + ih
                F_Ps_nloc[ipol,ia] -= 2.0*fac*Deff[ih,ih,ia] *
                   real( conj(dbetaNL_psi[ikb,ist]) * betaNL_psi[ikb,ist] )
            end
            
            #if psp.is_ultrasoft  # || upf(nt)%is_multiproj # need is_multiproj?
            # this case is almost always true for our case
                for ih in 1:nh[isp]
                    ikb = ijkb0 + ih
                    # in US case there is a contribution for jh /= ih. 
                    # We use here the symmetry in the interchange  of ih and jh
                    for jh in (ih+1):nh[isp]
                        jkb = ijkb0 + jh
                        F_Ps_nloc[ipol,ia] -= 2.0*fac*Deff[ih,jh,ia] * 
                            real( conj(dbetaNL_psi[ikb,ist]) * betaNL_psi[jkb,ist] +
                                  dbetaNL_psi[jkb,ist] * conj(betaNL_psi[ikb,ist]) )
                    end
                end
            #end # is_ultrasoft
        end
    end
    return
end


# This routine computes the contribution to atomic forces due
# to the dependence of the Q function on the atomic position.
# \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
#    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
# where:
# \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
#    w_i\langle \beta_m|\psi_i\rangle \]
# On output: the contribution is added to \(\text{forcenl}\).
function _add_F_uspp!(atoms, pw::PWGrid, pspots, pspotNL, potentials, F_Ps_nloc)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    atpos = atoms.positions

    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)
    G = pw.gvec.G
    G2 = pw.gvec.G2

    lmaxkb = pspotNL.lmaxkb
    nh = pspotNL.nh
    becsum = pspotNL.becsum

    Nspin = size(potentials.Total, 2)

    ok_uspp_or_paw = any(pspotNL.are_ultrasoft) || any(pspotNL.are_paw)
    if !ok_uspp_or_paw
        return
    end

    planfw = plan_fft!( zeros(ComplexF64,pw.Ns) ) # using default plan

    F_uspp = zeros(Float64, 3, Natoms)
    fact = 1.0*pw.CellVolume # this is 2*pw.CellVolume if using gamma only

    # Fourier transform of the total effective potential
    Vg = zeros(ComplexF64, Ng, Nspin)
    aux = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        aux[:] .= potentials.Total[:,ispin]
        #R_to_G!(pw, aux) #XXX why this will segfault?
        ff = reshape(aux, pw.Ns)
        planfw*ff
        @views aux[:] /= Npoints # rescale
        # Note the factors -i and 2pi/a *units of G) here in V(G)
        for ig in 1:Ng
            ip = idx_g2r[ig]
            Vg[ig,ispin] = -im*aux[ip] # XXX need factor tpiba?
        end
    end
    # Finished calculation -im*V_eff(G)


    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    Ylm_real_qe!(_lmax, G, ylmk0) # Ylm_real_qe accept l value starting from 0
    
    for isp in 1:Nspecies

        psp = pspots[isp]
        need_augmentation = psp.is_ultrasoft || psp.is_paw
        if !need_augmentation
            continue # skip for this species
        end

        # nij = max number of (ih,jh) pairs per atom type nt
        nij = Int64(nh[isp]*(nh[isp] + 1)/2)
        Qgm = zeros(ComplexF64, Ng, nij)
        ijh = 0
        for ih in 1:nh[isp], jh in ih:nh[isp]
            ijh = ijh + 1
            @views qvan2!( pspotNL, ih, jh, isp, G2, ylmk0, Qgm[:,ijh] )
            #println("ijh = ", ijh, " sum(Qgm) = ", sum(Qgm[:,ijh]))
        end
        #
        # nab = number of atoms of type nt
        nab = sum(atm2species .== isp)
        #
        aux1 = zeros(ComplexF64, Ng, nab, 3)
        dDeeq = zeros(Float64, nij, nab, 3, Nspin)
        #
        for ispin in 1:Nspin
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                nb = nb + 1
                # aux1 = product of potential, structure factor and iG
                for ig in 1:Ng
                    GX = G[1,ig]*atpos[1,ia] + G[2,ig]*atpos[2,ia] + G[3,ig]*atpos[3,ia]
                    Sf = cos(GX) + im*sin(GX)
                    cfac = Vg[ig,ispin] * Sf
                    aux1[ig,nb,1] = G[1,ig] * cfac
                    aux1[ig,nb,2] = G[2,ig] * cfac
                    aux1[ig,nb,3] = G[3,ig] * cfac
                end
            end
            println("\nsum aux1 = ", sum(aux1))
            println("sum Qgm = ", sum(Qgm))
            # dDeeq = dot product of aux1 with the Q functions
            # No need for special treatment of the G=0 term (is zero)
            for ipol in 1:3
                #CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, &
                #    aux1(1,1,ipol), 2*ngm, 0.0_dp, ddeeq(1,1,ipol,is), nij )
                dDeeq[:,:,ipol,ispin] .= fact * real(Qgm' * aux1[:,:,ipol]) # XXX
            end
        end
        #
        for ispin in 1:Nspin
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue # skip
                end
                nb = nb + 1
                for ipol in 1:3, ijh in 1:nij
                    F_uspp[ipol,ia] += dDeeq[ijh,nb,ipol,ispin] * becsum[ijh,ia,ispin]
                end
            end
        end
    end
    # Add F_uspp to the output array
    F_Ps_nloc[:,:] .+= F_uspp[:,:]
    return
end
