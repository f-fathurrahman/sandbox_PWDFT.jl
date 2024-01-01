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



# This routine computes the contribution to atomic forces due
# to the dependence of the Q function on the atomic position.
# \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
#    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
# where:
# \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
#    w_i\langle \beta_m|\psi_i\rangle \]
# On output: the contribution is added to \(\text{forcenl}\).


function _add_F_uspp!(atoms, pspotNL, F_Ps_nloc)
    # REAL(DP) :: fact
    # COMPLEX(DP) :: cfac
    # ! work space
    # COMPLEX(DP), ALLOCATABLE :: aux(:), aux1(:,:,:), vg(:,:), qgm(:,:)
    # REAL(DP),    ALLOCATABLE :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:), forceq(:,:)

    #
    ok_uspp_or_paw = any(pspotNL.are_ultrasoft) || any(pspotNL.are_paw)

    if !ok_uspp_or_paw
        return
    end

    F_uspp = zeros(Float64, 3, Natoms)

    fact = 1.0 # this is 2 if using gamma only

    #
    # fourier transform of the total effective potential
    #
    vg = zeros(ComplexF64, Ng, Nspin)
    aux = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        aux[:] .= potentials.Total[:,ispin]
        R_to_G!(pw, aux)
        # Note the factors -i and 2pi/a *units of G) here in V(G)
        for ig in 1:Ng
            ip = idx_g2r[ig]
            vg[:,ispin] = -im*aux[ip] # XXX need factor tpiba?
        end
    end
    # Finish calculation -im*V_eff(G)


    ylmk0 = zeros(Float64, lmaxq*lmaxq)
    ALLOCATE( ylmk0(ngm, lmaxq*lmaxq) )
    CALL ylmr2( lmaxq * lmaxq, ngm, g[:,1:ngm), gg(1:ngm), ylmk0 )
    
    ALLOCATE( qmod(ngm) )
    for ig in 1:ngm
    qmod(ig) = SQRT( gg(ig) )*tpiba
  ENDDO
  !
  DO nt = 1, ntyp
    IF ( upf(nt)%tvanp ) THEN
      !
      ! nij = max number of (ih,jh) pairs per atom type nt
      ! qgm contains the Q functions in G space
      !
      nij = nh(nt)*(nh(nt)+1)/2
      ALLOCATE( qgm(ngm,nij) )
      ijh = 0
      DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
          ijh = ijh + 1
          CALL qvan2( ngm, ih, jh, nt, qmod, qgm(:,ijh), ylmk0 )
        ENDDO
      ENDDO
      !
      ! nab = number of atoms of type nt
      !
      nab = 0
      DO na = 1, nat
        IF( ityp(na) == nt ) nab = nab + 1
      ENDDO
      ALLOCATE( aux1(ngm, nab, 3) )
      ALLOCATE( ddeeq(nij, nab, 3, nspin_mag) )
      !
      DO is = 1, nspin_mag
        nb = 0
        DO na = 1, nat
          IF (ityp(na) == nt) THEN
            nb = nb + 1
            !
            ! aux1 = product of potential, structure factor and iG
            !
            DO ig = 1, ngm
              cfac = vg(ig, is) * &
                   CONJG(eigts1(mill(1,ig),na) * eigts2(mill(2,ig),na) * eigts3(mill(3,ig),na) )
              aux1(ig,nb,1) = g(1,ig) * cfac
              aux1(ig,nb,2) = g(2,ig) * cfac
              aux1(ig,nb,3) = g(3,ig) * cfac
            ENDDO
            !
          ENDIF
        ENDDO
        !
        ! ddeeq = dot product of aux1 with the Q functions
        ! No need for special treatment of the G=0 term (is zero)
        !
        DO ipol = 1, 3
          CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, &
               aux1(1,1,ipol), 2*ngm, 0.0_dp, ddeeq(1,1,ipol,is), nij )
        ENDDO
        !
      ENDDO
      !
      DEALLOCATE(aux1)
      DEALLOCATE(qgm)
      !
      DO is = 1, nspin_mag
        nb = 0
        DO na = 1, nat
          IF( ityp(na) == nt ) THEN
            nb = nb + 1
            DO ipol = 1, 3
              DO ijh = 1, nij
                 forceq(ipol,na) = forceq(ipol,na) + ddeeq(ijh,nb,ipol,is) * becsum(ijh,na,is)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE( ddeeq )
    ENDIF
  ENDDO
  !
  10 CONTINUE


  F_nl[:,:] .+= F_uspp[:,:]

  return

end