function calc_stress_Ps_nloc!( atoms, pw, pspots, pspotNL, electrons, psiks, stress_Ps_nloc )

    G = pw.gvec.G
    idx_gw2g = pw.gvecw.idx_gw2g
    Ngw = pw.gvecw.Ngw
    Ngwx = maximum(Ngw)
    kfac = ones(Float64, Ngwx) # assume qcutz=0
    Gk = zeros(Float64, Ngwx, 3)
    Gk_length_inv = zeros(Float64, Ngwx)
    #
    k = pw.gvecw.kpoints.k
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    #
    Nstates = electrons.Nstates
    Nspin = electrons.Nspin
    Focc = electrons.Focc

    fill!(stress_Ps_nloc, 0.0)
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin - 1)*Nkpt
        #
        for igw in 1:Ngw[ik]
            ig = idx_gw2g[ik][igw]
            for i in 1:3
                Gk[igw,i] = G[i,ig] + k[i,ik]
            end
            # NOTE: reversed the usual dimension of Gk[1:3,1:Ngwx] to Gk[1:Ngwx,1:3]
            q = sqrt( Gk[igw,1]^2 + Gk[igw,2]^2 + Gk[igw,3]^2 )
            if q > 1e-8
                Gk_length_inv[igw] = 1.0 / q
            else
                Gk_length_inv[igw] = 0.0
            end
        end

        _stress_Ps_nloc_k!(
            ispin, ik,
            atoms, pw, pspots, pspotNL, electrons,
            Gk, Gk_length_inv,
            psiks, stress_Ps_nloc)

    end


    for l in 1:3, m in 1:(l-1)
        stress_Ps_nloc[m,l] = stress_Ps_nloc[l,m]
    end
    # Scale by -1/pw.CellVolume
    stress_Ps_nloc[:,:] .*= (-1/pw.CellVolume)


    return
end


function _stress_Ps_nloc_k!(
    ispin::Int64, ik::Int64,
    atoms, pw, pspots, pspotNL, electrons, 
    Gk, Gk_length_inv,
    psiks, stress_Ps_nloc
)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    Nstates = electrons.Nstates
    ebands = electrons.ebands

    nh = pspotNL.nh
    nhm = pspotNL.nhm
    Deeq = pspotNL.Deeq
    indv_ijkb0 = pspotNL.indv_ijkb0
    NbetaNL = pspotNL.NbetaNL

    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    Ngwk = pw.gvecw.Ngw[ik]

    ikspin = ik + (ispin - 1)*Nkpt
    psi = psiks[ikspin]
    Vnl_KB = pspotNL.betaNL[ik]
    betaNL_psi = Vnl_KB' * psi  # XXX precalculate this

    fac = wk[ik] # XXX: this is different from QE, it does not depend on ist
    # XXX Probably Focc also should be accounted for (also some normalization by Nstates)

    Deff = zeros(Float64, nhm, nhm, Natoms)

    evps = 0.0
    for ist in 1:Nstates
        if abs(fac) < 1e-9
            continue
        end
        _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ikspin], Deff )
        #
        ijkb0 = 0
        for isp in 1:Nspecies, ia in 1:Natoms
            #
            if atm2species[ia] != isp
                continue
            end
            #
            psp = pspots[isp]
            #
            for ih in 1:nh[isp]
                ikb = ijkb0 + ih
                evps += fac * Deff[ih,ih,ia] * abs(betaNL_psi[ikb,ist])^2
                #
                if psp.is_ultrasoft || (psp.Nproj > 1)
                    # only in the US case there is a contribution for jh != ih
                    # we use here the symmetry in the interchange of ih and jh
                    for jh in (ih + 1):nh[isp]
                        jkb = ijkb0 + jh
                        bibj = real( conj(betaNL_psi[ikb,ist]) * betaNL_psi[jkb,ist] )
                        evps += Deff[ih,jh,ia] * fac * bibj
                    end
                end
            end
            ijkb0 += nh[isp]
        end
    end # Nstates
    #    
    for l in 1:3
        stress_Ps_nloc[l,l] -= evps
    end


    #
    # non diagonal contribution - derivative of the bessel function
    #
    dvkb = zeros(ComplexF64, Ngwk, NbetaNL)
  
    _gen_us_dj!(ik, atoms, pw, pspots, pspotNL, dvkb)
    work2 = zeros(ComplexF64, Ngwk)
    work1 = zeros(ComplexF64, Ngwk)
    #
    for ist in 1:Nstates
        #
        fill!(work2, 0.0 + im*0.0)
        #
        _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ikspin], Deff )
        #
        ijkb0 = 0
        for isp in 1:Nspecies, ia in 1:Natoms
            #
            if atm2species[ia] != isp
                continue
            end
            #
            psp = pspots[isp]
            uspp_or_multiproj = psp.is_ultrasoft || (psp.Nproj > 1)
            #
            for ih in 1:nh[isp]
                ikb = ijkb0 + ih
                if !uspp_or_multiproj
                    ps = betaNL_psi[ikb,ist] * Deeq[ih,ih,ia,ispin]
                else
                    # in the USPP case there is a contribution also for jh != ih
                    ps = 0.0 + im*0.0
                    for jh in 1:nh[isp]
                        jkb = ijkb0 + jh
                        ps += betaNL_psi[jkb,ist] * Deff[ih,jh,ia]
                    end
                end
                #
                for igw in 1:Ngwk
                    work2[igw] = ps * dvkb[igw,ikb] + work2[igw]
                end
            end
            ijkb0 = ijkb0 + nh[isp]
        end
        #
        for ipol in 1:3, jpol in 1:ipol
            for igw in 1:Ngwk
                work1[igw] = psi[igw,ist] * Gk[igw,ipol] * Gk[igw,jpol] * Gk_length_inv[igw]
            end
            dd = real(dot(work1, work2))
            stress_Ps_nloc[ipol,jpol] -= wk[ik] * dd
        end
    end # ist


    #
    # non diagonal contribution - derivative of the spherical harmonics
    # (no contribution from l=0)

    # Immediate return
    if pspotNL.lmaxkb == 0
        return
    end

    xyz = zeros(Float64, 3, 3)
    xyz[1,:] = [1.0, 0.0, 0.0]
    xyz[2,:] = [0.0, 1.0, 0.0]
    xyz[3,:] = [0.0, 0.0, 1.0]
    # also can use Matrix(I(3))

    for ipol in 1:3

        _gen_us_dy!(ik, xyz[:,ipol], atoms, pw, pspots, pspotNL, dvkb)
        
        for ist in  1:Nstates
            #
            fill!(work2, 0.0 + im*0.0)
            _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ikspin], Deff )

            ijkb0 = 0
            for isp in 1:Nspecies, ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                psp = pspots[isp]
                uspp_or_multiproj = psp.is_ultrasoft || (psp.Nproj > 1)
                for ih in 1:nh[isp]
                    ikb = ijkb0 + ih
                    if !uspp_or_multiproj
                        ps = betaNL_psi[ikb,ist] * Deeq[ih,ih,ia,ispin]
                    else
                        # in the US case there is a contribution also for jh<>ih
                        ps = 0.0 + im*0.0
                        for jh in 1:nh[isp]
                            jkb = ijkb0 + jh
                            ps += betaNL_psi[jkb,ist] * Deff[ih,jh,ia]
                        end
                    end
                    #
                    for igw in 1:Ngwk
                        work2[igw] = ps * dvkb[igw,ikb] + work2[igw]
                    end
                end
                ijkb0 += nh[isp]
            end
            for jpol in 1:ipol
                for igw in 1:Ngwk
                    work1[igw] = psi[igw,ist] * Gk[igw,jpol]
                end
                dd = real(dot(work1, work2))
                stress_Ps_nloc[ipol,jpol] -= wk[ik] * dd
            end
        end # ist
    end # ipol
    #
    return

end