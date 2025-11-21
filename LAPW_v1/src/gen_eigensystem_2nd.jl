
# from eveqnsv
function gen_eigensystem_2nd!(
    ik::Int64,
    atoms::Atoms, pw::PWGrid, mt_vars::MuffinTins, apwlo_vars,
    apwalm, elec_chgst, bsmt, bsir, ndmag
)
    # XXX HARDCODED!!!
    ncmag = false
    spinorb = false

    nspinor = elec_chgst.nspinor
    spinpol = elec_chgst.spinpol
    nstsv = elec_chgst.nstsv
    nstfv = elec_chgst.nstfv

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    npcmt = mt_vars.npcmt
    npcmtmax = maximum(npcmt)
    #nrcmt = mt_vars.nrcmt
    #nrcmti = mt_vars.nrcmti
    #npcmti = mt_vars.npcmti
    lradstp = mt_vars.lradstp
    wrcmt = mt_vars.wrcmt

    # Read 1st variational eigenvalues and eigenvectors from files
    # XXX: Read this or just make them into arguments ?
    evalfv = deserialize("evals_1st_ik_$(ik).dat")
    evecfv = deserialize("evecs_1st_ik_$(ik).dat")

    evecsv = zeros(ComplexF64, nstsv, nstsv)
    @views evalsvp = elec_chgst.evalsv[:,ik]

    tevecsv = false
    if spinpol
        tevecsv = true
    end
    # TODO: Hartree-Fock/RDMFT/TDDFT/GW requires second-variational eigenvectors

    #
    # no calculation of second-variational eigenvectors
    if !tevecsv
        # evalsvp = evalfv
        for i in 1:nstsv
            evalsvp[i] = evalfv[i]
        end
        #
        # evecsv is identity matrix (ComplexF64)
        fill!(evecsv, 0.0)
        for i in 1:nstsv
            evecsv[i,i] = 1.0
        end
        return
    end
    #
    # number of spin combinations after application of Hamiltonian
    if spinpol
        if ncmag || spinorb
            nsc = 3
        else
            nsc = 2
        end
        nsd = 2
    else
        # XXX ffr: In what case we reach here?
        nsc = 1
        nsd = 1
    end
    #
    # special case of spin-orbit coupling and collinear magnetism
    if spinorb && cmagz
        socz = true
    else
        socz = false
    end
    #println("nsc = $nsc nsd = $nsd socz = $socz")
    #
    # zero the second-variational Hamiltonian (stored in the eigenvector array)
    fill!( evecsv, 0.0 )

    #
    # -------- muffin-tin part ---------------
    #
    wfmt1 = zeros(ComplexF64, npcmtmax, nstfv)
    wfmt2 = zeros(ComplexF64, npcmtmax)
    wfmt3 = zeros(ComplexF64, npcmtmax)
    wfmt4 = zeros(ComplexF64, npcmtmax, nsc)
    #
    #if (spinorb) allocate(zlflm(lmmaxo,3))
    #if (tafield .or. (xcgrad == 4)) allocate(gzfmt(npcmtmax,3))
    #
    # begin loop over atoms
    for ia in 1:Natoms
        isp = atm2species[ia]
        #nrc = nrcmt[isp]
        #nrci = nrcmti[isp]
        #nrco = nrc - nrci
        #irco = nrci + 1
        npc = npcmt[isp]
        #npci = npcmti[isp]
        #
        # compute the first-variational wavefunctions
        for ist in 1:nstfv
            #println("sum evecfv[:,ist] = ", sum(evecfv[:,ist]))
            @views wavefmt!( lradstp, ia, atoms, mt_vars, apwlo_vars,
                pw.gvecw.Ngw[ik], apwalm[ia], evecfv[:,ist], wfmt1[:,ist] )
            #println("sum wfmt1[:,ist] = ", sum(wfmt1[:,ist]))
        end
        #println("sum(wfmt1) = ", sum(wfmt1))
        #println("sum bsmt[ia][1:npc,ndmag] = ", sum(bsmt[ia][1:npc,ndmag]))
        #
        # begin loop over states
        for jst in 1:nstfv
            if spinpol
                # convert wavefunction to spherical coordinates
                @views backward_SHT!(mt_vars, isp, wfmt1[:,jst], wfmt2, coarse=true)
                #
                # apply Kohn-Sham effective magnetic field
                wfmt3[1:npc] .= bsmt[ia][1:npc,ndmag] .* wfmt2[1:npc]  # ffr: THIS IS IMPORTANT
                #
                # convert to spherical harmonics and store in wfmt4
                @views forward_SHT!(mt_vars, isp, wfmt3, wfmt4[:,1], coarse=true)
                wfmt4[1:npc,2] .= -wfmt4[1:npc,1]
                #
                # non-collinear magnetic field
                #if (ncmag) then
                #    !write(*,*) 'ncmag is applied'
                #    wfmt3(1:npc) = cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc)
                #    call zfsht(nrc, nrci, wfmt3, wfmt4(:,3))
                #endif
                if socz
                    @. wfmt4[1:npc,3] = 0.0
                end
                #
                # apply spin-orbit coupling if required
                # SKIPPED
            else
                # No SOC, no spinpol (HF case??)
                for k in 1:nsc
                    wfmt4[1:npc,k] = 0.0
                end
            end # if
            #
            # apply muffin-tin potential matrix if required
            # SKIPPED
            #
            # apply vector potential if required
            # SKIPPED
            #
            # add to second-variational Hamiltonian matrix
            #
            # upper diagonal block
            for ist in 1:jst
                z1 = zfmtinp(mt_vars, isp, wrcmt[isp], wfmt1[:,ist], wfmt4[:,1], coarse=true)
                #println("Upper diagonal block: ist=$ist z1 from zfmtinp = $z1")
                evecsv[ist,jst] += z1
            end
            #
            # lower diagonal block
            if nsc >= 2 # spinpol or ncmag or SOC
                j = jst + nstfv
                for ist in 1:jst
                    i = ist + nstfv
                    z1 = zfmtinp(mt_vars, isp, wrcmt[isp], wfmt1[:,ist], wfmt4[:,2], coarse=true)
                    evecsv[i,j] += z1
                end
            end
            #
            # off-diagonal block
            if nsc == 3 # in case of ncmag .or. spinorb
                @error "Should not go here"
                for ist in 1:nstfv
                    z1 = zfmtinp(mt_vars, isp, wrcmt[:,isp], wfmt1[:,ist], wfmt4[:,3], coarse=true)
                    evecsv[ist,j] += z1
                end
            end # if
        end # over states
        #
        # apply tau-DFT non-multiplicative potential if required
        # SKIPPED
    end # end loop over atoms
    
    #deallocate(wfmt2,wfmt3,wfmt4)
    #if (spinorb) deallocate(zlflm)
    #if (tafield.or.(xcgrad == 4)) deallocate(gzfmt)
    #deallocate(wfmt1)
    #if (xcgrad == 4) deallocate(gwfmt)
  
    idx_g2r = pw.gvec.idx_g2r
    idx_gw2g = pw.gvecw.idx_gw2g
    Ngw = pw.gvecw.Ngw
    Npoints = prod(pw.Ns)
    #
    # interstitial part     !
    if spinpol # || tafield || (xcgrad == 4)) then
        if socz
            nsc = 2
        end        
        wfir1 = zeros(ComplexF64, Npoints)
        wfir2 = zeros(ComplexF64, Npoints)
        wfgp = zeros(ComplexF64, Ngw[ik], nsc)
        # begin loop over states
        for jst in 1:nstfv
            fill!( wfir1, 0.0)
            for igw in 1:Ngw[ik]
                ig = idx_gw2g[ik][igw]
                ip = idx_g2r[ig]
                wfir1[ip] = evecfv[igw,jst]
            end
            #
            # Fourier transform wavefunction to real-space
            G_to_R!(pw, wfir1)
            #
            # multiply with magnetic field and transform to G-space
            if spinpol
                wfir2[:] = bsir[:,ndmag] .* wfir1 # IMPORTANT: multiply in real-space
                R_to_G!(pw, wfir2)
                for igw in 1:Ngw[ik]
                    ig = idx_gw2g[ik][igw]
                    ip = idx_g2r[ig]
                    wfgp[igw,1] = wfir2[ip]
                end
                wfgp[1:Ngw[ik],2] = -wfgp[1:Ngw[ik],1]
                #
                #if ncmag
                #    wfir2(:) = cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
                #    #call zfftifc(3,ngridg,-1,wfir2)
                #    for igp in 1:ngp
                #        wfgp(igp,3)=wfir2(igfft(igpig(igp)))
                #    end
                #end
            else
                wfgp[1:Ngw[ik],1:nsd] = 0.0
            end
            #
            # apply vector potential if required
            # SKIPPED
            #
            # add to second-variational Hamiltonian matrix
            #
            # upper diagonal block
            for ist in 1:jst
                evecsv[ist,jst] += BLAS.dotc(Ngw[ik], evecfv[:,ist], 1, wfgp[:,1], 1)
            end
            #
            # lower diagonal block
            if nsc >= 2
                j = jst + nstfv
                for ist in 1:jst
                    i = ist + nstfv
                    evecsv[i,j] += BLAS.dotc(Ngw[ik], evecfv[:,ist], 1, wfgp[:,2], 1)
                end
            end
            #
            # off-diagonal block
            if nsc == 3
                for ist in 1:nstfv
                    evecsv[ist,j] += BLAS.dotc(Ngw[ik], evecfv[:,ist], 1, wfgp[:,3], 1)
                    # set conjugate?
                end
            end
        end    # end loop over states nstfv
        #deallocate(wfir1, wfir2, wfgp)
    end
    #
    # add the diagonal first-variational part
    i = 0
    #println("nspinor = ", nspinor)
    for ispn in 1:nspinor
        for ist in 1:nstfv
            i += 1
            evecsv[i,i] += evalfv[ist]
        end
    end
    #
    #println("sum evecsv before diagonalization = ", sum(evecsv))
    serialize("Ham_2nd_ik_$(ik).dat", evecsv)
    #
    spcpl = false
    if spcpl || (!spinpol)
        # spins are coupled; or spin-unpolarised: full diagonalisation
        evalsvp[:], evecsv[:,:] = eigen(Hermitian(evecsv))
    else
        # spins not coupled: block diagonalize H
        idx = 1:nstfv
        evalsvp[idx], evecsv[idx,idx] = eigen(Hermitian(evecsv[idx,idx]))
        i = nstfv + 1
        idx = (nstfv+1):nstsv
        evalsvp[idx], evecsv[idx,idx] = eigen(Hermitian(evecsv[idx,idx]))
        for i in 1:nstfv, j in 1:nstfv
            evecsv[i, j+nstfv] = 0.0
            evecsv[i+nstfv, j] = 0.0
        end
    end

    serialize("evals_2nd_ik_$(ik).dat", evalsvp)
    serialize("evecs_2nd_ik_$(ik).dat", evecsv)

    return

end
