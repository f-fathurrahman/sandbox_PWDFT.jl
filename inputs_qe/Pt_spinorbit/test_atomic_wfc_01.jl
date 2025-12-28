function calc_Natomwfc( atoms, pspots, noncollinear )
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Natomwfc = 0
    SMALL_OCC = 1e-12
    for ia in 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for n in 1:psp.Nchi
            if psp.occ_chi[n] < SMALL_OCC
                continue
            end
            if noncollinear
                if psp.has_so
                    Natomwfc += 2*psp.lchi[n]
                    if abs(psp.jchi[n] - psp.lchi[n] - 0.5) < 1e-6
                        Natomwfc += 2
                    end
                else
                    Natomwfc += 2*( 2 * psp.lchi[n] + 1 )
                end # if
            else
                # The usual case (no spinorb and no noncollinear magn)
                Natomwfc += 2*psp.lchi[n] + 1
            end
        end
    end
    return Natomwfc
end




function debug_atomic_wfc(Ham::Hamiltonian)
    atoms = Ham.atoms
    pspots = Ham.pspots
    pw = Ham.pw
    noncollinear = Ham.electrons.noncollinear
    #
    ik = 1
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atpos = atoms.positions
    atm2species = atoms.atm2species
    Ngwk = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k
    rot_ylm = Ham.pspotNL.rot_ylm
    # calculate max angular momentum required in wavefunctions
    lchi_max = 0
    Nchi_max = 0
    for isp in 1:Nspecies
        psp = pspots[isp]
        lchi_max = max(lchi_max, maximum(psp.lchi[1:psp.Nchi]))
        Nchi_max = max(Nchi_max, psp.Nchi)
    end

    Gk = zeros(Float64, 3, Ngwk)
    Gk_len = zeros(Float64, Ngwk)
    # 
    # igk: index for Gk
    # ig: index for G
    for igk in 1:Ngwk
        ig = idx_gw2g[ik][igk] # index of Gk in G
        Gk[1,igk] = G[1,ig] + k[1,ik]
        Gk[2,igk] = G[2,ig] + k[2,ik]
        Gk[3,igk] = G[3,ig] + k[3,ik]
        Gk_len[igk] = sqrt(Gk[1,igk]^2 +  Gk[2,igk]^2 + Gk[3,igk]^2)
    end

    ylm = zeros(Float64, Ngwk, (lchi_max+1)^2)
    # Ylm_real_qe accept l value starting from 0 (the actual 'physics' angular momentum number)
    Ylm_real_qe!(lchi_max, Gk, ylm)

    dq = 0.01 # XXX HARDCODED
    n_starting_wfc = 0
    # chiq = radial fourier transform of atomic orbitals chi
    SMALL_OCC = 1e-12
    chiq = zeros(Float64, Ngwk, Nchi_max, Nspecies)
    for isp in 1:Nspecies
        psp = pspots[isp]
        tab_at = PWDFT._build_chi_interp_table(psp, pw)
        for nb in 1:psp.Nchi
            # skip if occ_chi is very small (it is originally zero in pw.x)
            if psp.occ_chi[nb] < SMALL_OCC
                continue
            end
            for igk in 1:Ngwk
                px = Gk_len[igk]/dq - floor(Int64, Gk_len[igk]/dq)
                ux = 1.0 - px
                vx = 2.0 - px
                wx = 3.0 - px
                i0 = floor(Int64, Gk_len[igk] / dq ) + 1
                i1 = i0 + 1
                i2 = i0 + 2
                i3 = i0 + 3
                chiq[igk,nb,isp] = tab_at[i0,nb] * ux * vx * wx / 6.0 +
                                   tab_at[i1,nb] * px * vx * wx / 2.0 -
                                   tab_at[i2,nb] * px * ux * wx / 2.0 +
                                   tab_at[i3,nb] * px * ux * vx / 6.0
            end
        end
    end

    Natomwfc = calc_Natomwfc(atoms, pspots, noncollinear)

    Ngwx = maximum(pw.gvecw.Ngw)
    Npol = 1
    if noncollinear
        Npol = 2
    end
    n_starting_wfc = max(Natomwfc, Ham.electrons.Nstates )
    wfcatom = zeros(ComplexF64, Ngwx, Npol, n_starting_wfc )

    aux = zeros(ComplexF64, Ngwk)
    Sf = zeros(ComplexF64, Ngwk)
    fill!(wfcatom, 0.0 + 0.0)
    n_starting_wfc = 0  # this should be modified

    starting_spin_angle = true # FIXME

    for ia in 1:Natoms
        #    
        for igk in 1:Ngwk
            # XXX use dot product here?
            GkX = atpos[1,ia]*Gk[1,igk] + atpos[2,ia]*Gk[2,igk] + atpos[3,ia]*Gk[3,igk]
            Sf[igk] = cos(GkX) - im*sin(GkX)
        end
        isp = atm2species[ia]
        psp = pspots[isp]
        for nb in 1:psp.Nchi
            if psp.occ_chi[nb] < SMALL_OCC
                continue
            end
            l = psp.lchi[nb]
            lphase = im^l
            #  the factor i^l MUST BE PRESENT in order to produce
            #  wavefunctions for k=0 that are real in real space
            if noncollinear
                if psp.has_so
                    if starting_spin_angle # || .NOT. domag) THEN
                        # spin-orbit coupling only, no magnetism
                        #CALL my_atomic_wfc_so()
                        n_starting_wfc = _atomic_wfc_so!(
                            nb, isp, l, n_starting_wfc, psp, Natomwfc, Ngwk, lphase,
                            chiq, rot_ylm, ylm, Sf, aux, wfcatom
                        )
                    else
                        # spin-orbit coupling with magnetism
                        #CALL my_atomic_wfc_so_mag()
                        error("Not yet implemented 133")
                    end
                else
                    # noncollinear magn, no spinpol
                    #CALL my_atomic_wfc_nc()
                    error("Not yet implemented 138")
                end
            else
                #CALL my_atomic_wfc___()
                error("Not yet implemented 142")
            end
        end # loop over nwfc
    end # loop over natoms

    @infiltrate

end

function _atomic_wfc_so!(
    nb, isp, l, n_starting_wfc, psp, Natomwfc, Ngwk, lphase,
    chiq, rot_ylm, ylm, Sf, aux, wfcatom
)
 
    lmaxx  = 3  # HARDCODED !!!!!
    # XXX Get this into pspotNL ????

    j = psp.jchi[nb]
    fact = zeros(Float64, 2)
    for m in (-l-1):l    
        fact[1] = PWDFT.spinor_coef(l, j, m, 1)
        fact[2] = PWDFT.spinor_coef(l, j, m, 2)
        if abs(fact[1]) > 1e-8 || abs(fact[2]) > 1e-8
            n_starting_wfc += 1
            if n_starting_wfc > Natomwfc
                error("atomic_wfc_so internal error: too many wfcs")
            end
            for is in 1:2
                if abs(fact[is]) > 1e-8
                    ind = lmaxx + 1 + PWDFT.sph_ind(l, j, m, is)
                    fill!(aux, 0.0 + im*0.0)
                    for n1 in 1:(2*l+1)
                        ind1 = l^2 + n1
                        if abs(rot_ylm[ind,n1]) > 1e-8
                            aux[:] .+= rot_ylm[ind,n1] * ylm[:,ind1]
                        end
                    end
                    for igk in 1:Ngwk
                        wfcatom[igk,is,n_starting_wfc] = lphase*fact[is]*Sf[igk]*aux[igk]*chiq[igk,nb,isp]
                    end
                else
                    @views wfcatom[:,is,n_starting_wfc] .= 0.0 + im*0.0
                end
            end
        end
    end
    return n_starting_wfc
end


