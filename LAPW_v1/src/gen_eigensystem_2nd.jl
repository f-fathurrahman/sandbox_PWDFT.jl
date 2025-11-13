
# from eveqnsv
function gen_eigensystem_2nd!(ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)

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
            evalsv[i] = evalfv[i]
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
    # coupling constant of the external A-field (1/c)
    ca = 1.0/solsc
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
    println("nsc = $nsc nsd = $nsd socz = $socz")
    #
    #ld = lmmaxdm*nspinor # probably not yet needed
    #
    # zero the second-variational Hamiltonian (stored in the eigenvector array)
    fill!( evecsv, 0.0)

    #
    # -------- muffin-tin part ---------------
    #
    wfmt1 = zeros(ComplexF64, npcmtmax,nstfv)
    #
    #if (xcgrad == 4) allocate(gwfmt(npcmtmax,3,nstfv))
    #
    wfmt2 = zeros(ComplexF64, npcmtmax)
    wfmt3 = zeros(ComplexF64, npcmtmax)
    wfmt4 = zeros(ComplexF64, npcmtmax,nsc)
    #
    #if (spinorb) allocate(zlflm(lmmaxo,3))
    #if (tafield .or. (xcgrad == 4)) allocate(gzfmt(npcmtmax,3))
    #
    # begin loop over atoms
    for ia in 1:Natoms
        isp = atm2species[ia]
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        nrco = nrc - nrci
        irco = nrci + 1
        npc = npcmt[isp]
        npci = npcmti[isp]
        #
        # compute the first-variational wavefunctions
        for ist in 1:nstfv
            wavefmt!( lradstp, ias, ngp, apwalm[ia], evecfv[:,ist], wfmt1[:,ist] )
        end
        #
        # begin loop over states
        for jst in 1:nstfv
            if spinpol
                # convert wavefunction to spherical coordinates
                backward_SHT(mt_vars, wfmt1[:,jst], wfmt2, coarse=true)
                #
                # apply Kohn-Sham effective magnetic field
                @. wfmt3[1:npc] = bsmt[ia][1:npc,ndmag]*wfmt2[1:npc]  # ffr: THIS IS IMPORTANT
                #
                # convert to spherical harmonics and store in wfmt4
                forward_SHT!(mt_vars, wfmt3, wfmt4, coarse=true)
                @. wfmt4[1:npc,2] = -wfmt4[1:npc,1]
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
                z1 = zfmtinp(mt_vars, wrcmt[:,isp], wfmt1[:,ist], wfmt4, coarse=true)
                evecsv[ist,jst] += z1
            end
            #
            # lower diagonal block
            if nsc >= 2 # spinpol or ncmag or SOC
                j = jst + nstfv
                for ist in 1:jst
                    i = ist + nstfv
                    z1 = zfmtinp(mt_vars, wrcmt[:,isp], wfmt1[:,ist], wfmt4[:,2], coarse=true)
                    evecsv[i,j] += z1
                end
            end
            #
            # off-diagonal block
            if nsc == 3 # in case of ncmag .or. spinorb
                for ist in 1:nstfv
                    z1 = zfmtinp(mt_vars, wrcmt[:,isp], wfmt1[:,ist], wfmt4[:,3], coarse=true)
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
  

    #     interstitial part     !
    if spinpol # || tafield || (xcgrad == 4)) then
        if socz
            nsc = 2
        end
        
        #if (xcgrad == 4) allocate(gwfgp(ngkmax,nstfv))
        
        wfir1 = zeros(ComplexF64, Npoints)
        wfir2 = zeros(ComplexF64, Npoints)
        wfgp = zeros(ComplexF64, ngkmax, nsc)
        # begin loop over states
        for jst in 1:nstfv
            fill!( wfir1, 0.0)
            for ig in 1:Ng
                wfir1(igfft(igpig(igp)))=evecfv(igp,jst)
            end
            #
            # Fourier transform wavefunction to real-space
            G_to_R!(pw, wfir1)
            #call zfftifc(3,ngridg,1,wfir1)
            #
            # multiply with magnetic field and transform to G-space
            if spinpol
                wfir2[:] = bsir[:,ndmag] .* wfir1 # IMPORTANT: multiply in real-space
                #call zfftifc(3, ngridg, -1, wfir2)
                R_to_G!(pw, wfir2)
                for igp in 1:ngp
                    wfgp[igp,1] = wfir2(igfft(igpig(igp)))
                end
                wfgp[1:ngp,2] = -wfgp[1:ngp,1]
                #
                if ncmag
                    wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
                    #call zfftifc(3,ngridg,-1,wfir2)
                    for igp in 1:ngp
                        wfgp(igp,3)=wfir2(igfft(igpig(igp)))
                    end
                end
            else
                wfgp[1:ngp,1:nsd] = 0.0
            end
            #
            # apply vector potential if required
            # SKIPPED
            #
            # add to second-variational Hamiltonian matrix
            #
            # upper diagonal block
            for ist in 1:jst
                evecsv[ist,jst] = BLAS.dotc(ngp, evecfv[:,ist], 1, wfgp, 1)
            end
            #
            # lower diagonal block
            if nsc >= 2
                j = jst + nstfv
                for ist in 1:jst
                    i = ist + nstfv
                    evecsv[i,j] += BLAS.dotc(ngp, evecfv[:,ist], 1, wfgp[:,2], 1)
                end
            end
            #
            # off-diagonal block
            if nsc == 3
                for ist in 1:nstfv
                    evecsv[ist,j] += BLAS.dotc(ngp, evecfv[:,ist], 1, wfgp[:,3], 1)
                    # set conjugate?
                end
            end
        end    # end loop over states nstfv
        #deallocate(wfir1, wfir2, wfgp)
    end
    #
    # add the diagonal first-variational part
    i = 0
    for ispn in 1:nspinor
        for ist in 1:nstfv
            i += 1
            evecsv[i,i] += evalfv[ist]
        end
    end
    #
    if spcpl || (!spinpol)
        # spins are coupled; or spin-unpolarised: full diagonalisation
        eveqnz(nstsv, nstsv, evecsv, evalsvp)
    else
        # spins not coupled: block diagonalise H
        eveqnz!( nstfv, nstsv, evecsv, evalsvp)
        i = nstfv + 1
        eveqnz!( nstfv, nstsv, evecsv[i,i], evalsvp[i])
        for i in 1:nstfv, j in 1:nstfv
            evecsv[i, j+nstfv] = 0.0
            evecsv[i+nstfv, j] = 0.0
        end
    end

    return

end
