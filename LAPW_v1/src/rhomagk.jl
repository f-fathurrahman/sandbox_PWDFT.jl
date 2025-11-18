
function rhomagk!(
    ik, atoms, pw, mt_vars, apwlo_vars, elec_chgst,
    apwalm, evecfv, evecsv,
    rhomt, rhoir;
    magmt=nothing, magir=nothing
)

    println("\n----- ENTER rhomagk for ik=$ik")

    epsocc = 1e-8

    # FIXME: need magir and magmt for spinpol
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    npcmt = mt_vars.npcmt
    npcmtmax = maximum(npcmt)
    lradstp = mt_vars.lradstp
    
    Ngw = pw.gvecw.Ngw
    wppt = pw.gvecw.kpoints.wk[ik]
    idx_gw2r = pw.gvecw.idx_gw2r
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume

    nstfv = elec_chgst.nstfv
    nstsv = elec_chgst.nstsv
    occsv = elec_chgst.occsv
    spinpol = elec_chgst.spinpol
    nspinor = elec_chgst.nspinor
    ncmag = false  # HARDCODED
    tevecsv = elec_chgst.tevecsv

    done = zeros(Bool, nstfv)
    if tevecsv
        wfmt1 = zeros(ComplexF64, npcmtmax, nstfv)
    end
    wfmt2 = zeros(ComplexF64, npcmtmax)
    wfmt3 = zeros(ComplexF64, npcmtmax, nspinor)

    # loop over all atoms
    for ia in 1:Natoms
        isp = atm2species[ia]
        npc = npcmt[isp]
        # de-phasing factor for spin-spirals
        # SKIPPED
        fill!(done, false)
        # loop over all second-variational state
        for j in 1:nstsv
            wo = occsv[j,ik] # occupation number
            # skip this state if it is empty or nearly empty
            if abs(wo) < epsocc
                continue
            end
            #
            wo = wo*wppt
            #
            # if using 2nd variational scheme
            if tevecsv
                # generate spinor wavefunction from second-variational eigenvectors
                i = 0
                for ispn in 1:nspinor
                    #jspn = jspnfv[ispn] # XXX jspn is always 1
                    wfmt3[1:npc,ispn] .= 0.0
                    for ist in 1:nstfv
                        i = i + 1
                        z1 = evecsv[i,j]
                        if abs(real(z1)) + abs(imag(z1)) > epsocc
                            if !done[ist]
                                @views wavefmt!(lradstp, ia, atoms, mt_vars, apwlo_vars, Ngw[ik], apwalm[ia], evecfv[:,ist], wfmt2)
                                # convert to spherical coordinates
                                @views backward_SHT!(mt_vars, isp, wfmt2, wfmt1[:,ist], coarse=true)
                                done[ist] = true
                            end # if
                            # add to spinor wavefunction
                            #CALL zaxpy(npc,z1,wfmt1(:,ist,jspn),1,wfmt3(:,ispn),1)
                            @. wfmt3[1:npc,ispn] += z1 * wfmt1[1:npc,ist]
                        end
                    end # do
                end # do
            else 
                # not using 2nd variational scheme
                #
                # spin-unpolarised wavefunction
                @views wavefmt!(lradstp, ia, atoms, mt_vars, apwlo_vars, Ngw[ik], apwalm[ia], evecfv[:,j], wfmt2)
                # The result is stored in wfmt2
                #
                # convert to spherical coordinates
                @views backward_SHT!(mt_vars, isp, wfmt2, wfmt3[:,1], coarse=true)
            end
            #
            # add to density and magnetisation
            if spinpol
                # spin-polarised
                if ncmag
                    # non-collinear
                    #rhomagk_rmk1(npc, wo, wfmt3, wfmt3[:,2], rhomt[:,ias], magmt[:,ias,1], magmt[:,ias,2], magmt[:,ias,3])
                    println("SHOULD NOT PASS HERE")
                else
                    # collinear
                    #println("sum rhomt before: ", sum(rhomt[ia]))
                    #println("sum wfmt3 = ", sum(wfmt3))
                    @views rhomagk_rmk2!(npc, wo, wfmt3[:,1], wfmt3[:,2], rhomt[ia], magmt[ia][:,1])
                    #println("sum rhomt after: ", sum(rhomt[ia]))
                end
            else
                # spin-unpolarized
                @views rhomagk_rmk3!(npc, wo, wfmt3[:,1], rhomt[ia])
            end 
  
        end # over states
  
    end # over atoms


    #------------------------------------------------!
    #     interstitial density and magnetisation     !
    #------------------------------------------------!
    #
    wfir = zeros(ComplexF64, Npoints, nspinor)
    #
    # loop over all states
    #
    for j in 1:nstsv
        wo = occsv[j,ik]
        if abs(wo) < epsocc
            #println("Skipped for j = $j")
            continue
        end
        wo = wo*wppt/CellVolume
        fill!(wfir, 0.0)
        #
        if tevecsv
            # generate spinor wavefunction from second-variational eigenvectors
            i = 0 # idx for accessing basis function in evecsv
            for ispn in 1:nspinor
                for ist in 1:nstfv
                    i = i + 1
                    z1 = evecsv[i,j]
                    if abs(real(z1)) + abs(imag(z1)) > epsocc # XXX why check with epsocc?
                        for igw in 1:Ngw[ik]
                            ip = idx_gw2r[ik][igw]
                            wfir[ip,ispn] += z1*evecfv[igw,ist]
                        end
                    end # if
                end
            end
        else
            # spin-unpolarised wavefunction
            for igw in 1:Ngw[ik]
                ip = idx_gw2r[ik][igw]
                wfir[ip,1] = evecfv[igw,j]
            end
        end
    
        # Fourier transform wavefunction to real-space
        #println("sum wfir before FFT: ", sum(wfir))
        for ispn in 1:nspinor
            @views G_to_R!(pw, wfir[:,ispn])
        end
        wfir *= Npoints # scale to match Elk convention
        #println("sum wfir after FFT: ", sum(wfir))
        # add to density and magnetisation
        if spinpol 
            # spin-polarised
            if ncmag
                # non-collinear
                @views rhomagk_rmk1!(Npoints, wo, wfir[:,1], wfir[:,2], rhoir, magir, magir[:,2], magir[:,3])
            else
                # collinear
                @views rhomagk_rmk2!(Npoints, wo, wfir[:,1], wfir[:,2], rhoir, magir)
            end
        else
            # XXX We pass full FFT grid array here, so ngtot -> Npoints
            # spin-unpolarised
            rhomagk_rmk3!(Npoints, wo, wfir, rhoir)
        end
    end # nstsv

    return
end


function rhomagk_rmk1!(n, wo, wf1, wf2, rho, mag1, mag2, mag3)
    for i in 1:n
        z1 = wf1[i]
        z2 = wf2[i]
        t1 = real(z1)^2 + imag(z1)^2
        t2 = real(z2)^2 + imag(z2)^2
        z1 = conj(z1)*z2
        mag1[i] = mag1[i] + wo2*dble(z1)
        mag2[i] = mag2[i] + wo2*aimag(z1)
        mag3[i] = mag3[i] + wo*(t1-t2)
        rho[i] += wo*(t1 + t2)
    end
    return
end


function rhomagk_rmk2!(n, wo, wf1, wf2, rho, mag)
    for i in 1:n
        t1 = real(wf1[i])^2 + imag(wf1[i])^2
        t2 = real(wf2[i])^2 + imag(wf2[i])^2
        mag[i] += wo*(t1 - t2)
        rho[i] += wo*(t1 + t2)
    end
    return
end

# used for non-spin-polarized case
function rhomagk_rmk3!(n, wo, wf, rho)
    for i in 1:n
        rho[i] += wo*( real(wf[i])^2 + imag(wf[i])^2 )
    end
    return
end