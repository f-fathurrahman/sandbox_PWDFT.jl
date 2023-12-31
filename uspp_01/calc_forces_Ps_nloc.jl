function my_calc_forces_Ps_nloc!(
    atoms::Atoms,    
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    electrons::Electrons,
    pspotNL::PsPotNL_UPF,
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
                dbetaNL[igw,ibeta] = -im * betaNL[ik][ig,ibeta] * G[ipol,ig]
            end
            # betaNL psi 
            @views dbetaNL_psi[:,:] .= dbetaNL[1:Ngw_ik,:]' * psi
            #
            # this will call sum over bands
            _force_Ps_nloc_k!(ipol, ik, ispin, atoms, pw, pspots,
                electrons, pspotNL, betaNL_psi, dbetaNL_psi, F_Ps_nloc)
        end
    end

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






function _calc_Deff!(
    ispin::Int64, atoms, pspotNL, ebnd::Float64, Deff
)
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    Deeq = pspotNL.Deeq
    qq_at = pspotNL.qq_at
    #
    ok_uspp_or_paw = any(pspotNL.are_ultrasoft) || any(pspotNL.are_paw)
    #
    @views Deff[:,:,:] .= Deeq[:,:,:,ispin]
    #println("ok_uspp_or_paw = ", ok_uspp_or_paw)
    #
    if ok_uspp_or_paw
        for isp in 1:Nspecies, ia in 1:Natoms
            if isp != atm2species[ia]
                continue
            end
            Deff[:,:,ia] .-= ebnd * qq_at[:,:,ia]
        end
    end
    #println("sum qq_at = ", sum(qq_at))
    #println("sum Deeq = ", sum(Deeq))
    println("sum Deff = ", sum(Deff))
    return
end