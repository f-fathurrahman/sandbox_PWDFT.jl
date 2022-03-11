function debug_rhoinit!(
    atoms, atsp_vars,
    mt_vars, pw,
    rhomt, rhoir
)
    # OUTPUT: rhomt and rhoir

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    ecutrho = pw.ecutrho
    gmaxvr = sqrt(2*ecutrho)
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    idx_g2r = pw.gvec.idx_g2r

    nrsp = atsp_vars.nrsp
    rsp = atsp_vars.rsp
    rhosp = atsp_vars.rhosp 

    # For debugging
    # println("Some rhosp: ")
    #for isp in 1:Nspecies
    #    println("rhosp for isp = ", isp)
    #    for ir in 1:10
    #        @printf("%4d %18.10e\n", ir, rhosp[isp][ir])
    #    end
    #end
    #exit()

    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    #
    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    #
    rcmt = mt_vars.rcmt
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    #
    npmt = mt_vars.npmt
    npcmt = mt_vars.npcmt

    # FIXME: need to precompute this?
    sfacg = calc_strfact(atoms, pw)

    # FIXME: need to precompute this?
    ylmg = zeros(ComplexF64, lmmaxo, Ng)
    genylmg!(mt_vars.lmaxo, pw.gvec.G, ylmg)

    lmax = min(mt_vars.lmaxi,1) # FIXME: Why using 1 instead of lmaxi

    epslat = 1e-6
  
    # zero the charge density arrays
    for ia in 1:Natoms
        fill!(rhomt[ia], 0.0)
    end
    fill!(rhoir, 0.0)
  
    #
    # compute the superposition of all the atomic density tails
    #   
    zfft = zeros(ComplexF64,Npoints)

    for isp in 1:Nspecies

        # local arrays inside the loop
        wr = zeros(Float64,nrsp[isp])
        fr = zeros(Float64,nrsp[isp])

        nr = nrmt[isp]
        nrs = nrsp[isp]
        nro = nrs - nr + 1

        # determine the weights for the radial integral
        # integrate for radial points outside the muffin-tin
        idx = nr:nr+nro-1 # or better: idx = nr:nrs
        @views wsplint!( nro, rsp[isp][idx] ,wr[idx] )
        #
        for ig in 1:Ng
            t1 = sqrt(G2[ig])
            # spherical bessel function j_0(x) times the atomic density tail
            if t1 > epslat
                # G != 0 term
                t2 = 1.0/t1
                for ir in nr:nrs
                  x = t1*rsp[isp][ir]
                  fr[ir] = t2*sin(x)*rhosp[isp][ir]*rsp[isp][ir]
                end
            else
                # G=0 term
                for ir in nr:nrs
                    fr[ir] = rhosp[isp][ir] * rsp[isp][ir]^2
                end
            end
            @views t1 = dot( wr[nr:nrs], fr[nr:nrs] )
            # apply low-pass filter
            t1 = t1*exp(-4.0*G2[ig]/gmaxvr^2)
            ffg = (4*pi/CellVolume)*t1
            #
            ip = pw.gvec.idx_g2r[ig]
            zfft[ip] = zfft[ip] + ffg*sfacg[ig,isp]  # do not use conj
        end
    end

    #println("Some zfft")
    #for ig in 1:10
    #    ip = pw.gvec.idx_g2r[ig]
    #    @printf("%8d %18.10e %18.10e\n", ip, real(zfft[ip]), imag(zfft[ip]))
    #end
    #exit()

    println("sum sfacg = ", sum(sfacg))
    s = 0.0
    for isp in 1:Nspecies
        s = s + sum(rhosp[isp])
    end
    println("sum rhosp = ", s)
    println("sum zfft  = ", sum(zfft))

    nrcmtmax = maximum(nrcmt)
    npcmtmax = maximum(npcmt)
    jl = OffsetArray(
        zeros(Float64,lmax+1,nrcmtmax), 0:lmax, 1:nrcmtmax
    )
    zfmt = zeros(ComplexF64, npcmtmax)
    
    for ia in 1:Natoms
        isp = atm2species[ia]
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        irco = nrci + 1
        zfmt[1:npcmt[isp]] .= 0.0
        println("nrci = ", nrci)
        println("nrc  = ", nrc)
        for ig in 1:Ng
            ip = idx_g2r[ig]
            for irc in 1:nrc
                x = sqrt(G2[ig])*rcmt[isp][irc]
                for l in 0:lmax
                    jl[l,irc] = sphericalbesselj(l, x)
                end
            end
            z1 = 4*pi*zfft[ip] * conj(sfacg[ig,isp]) # XXX using conj
            lm = 0
            for l in 0:lmax
                z2 = im^l * z1
                for m in -l:l
                    lm = lm + 1
                    z3 = z2*conj(ylmg[lm,ig])
                    i = lm
                    for irc in 1:nrci
                        zfmt[i] = zfmt[i] + jl[l,irc]*z3
                        i = i + lmmaxi
                    end
                    for irc in irco:nrc
                        zfmt[i] = zfmt[i] + jl[l,irc]*z3
                        i = i + lmmaxo
                    end
                end
            end
        end
        println("Finish loop over G")

        println("Some zfmt")
        for i in 1:10
            @printf("%8d %18.10e %18.10e\n", i, real(zfmt[i]), imag(zfmt[i]))
        end

        z_to_rf_mt!( mt_vars, nrc, nrci, zfmt, rhomt[ia] )

        println("Some rhomt after z_to_rf_mt")
        println("inner")
        for i in 1:10
            @printf("%8d %18.10e\n", i, rhomt[ia][i])
        end
        println("outer")
        for i in nrmti[isp]+1:nrmti[isp]+lmmaxo
            @printf("%8d %18.10e\n", i, rhomt[ia][i])
        end

    end

    # convert the density from a coarse to a fine radial mesh
    rf_mt_c_to_f!( atoms, atsp_vars, mt_vars, rhomt )

    isp = 1
    ia = 1
    println("Some rhomt after rf_mt_c_to_f")
    println("inner")
    for i in 1:10
        @printf("%8d %18.10e\n", i, rhomt[ia][i])
    end
    println("outer")
    for i in nrmti[isp]+1:nrmti[isp]+lmmaxo
        @printf("%8d %18.10e\n", i, rhomt[ia][i])
    end

    # add the atomic charge density and the excess charge in each muffin-tin
    chgexs = 0.0  # FIXME
    y00 = 0.28209479177387814347

    t1 = chgexs/CellVolume
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        i = 1
        for ir in 1:nri
            t2 = (t1 + rhosp[isp][ir])/y00
            rhomt[ia][i] = rhomt[ia][i] + t2
            i = i + lmmaxi
        end
        for ir in nri+1:nr
            t2 = (t1 + rhosp[isp][ir])/y00
            rhomt[ia][i] = rhomt[ia][i] + t2
            i = i + lmmaxo
        end
    end

    isp = 1
    ia = 1
    println("Some rhomt after adding rhosp")
    println("inner")
    for i in 1:10
        @printf("%8d %18.10e\n", i, rhomt[ia][i])
    end
    println("outer")
    for i in nrmti[isp]+1:nrmti[isp]+lmmaxo
        @printf("%8d %18.10e\n", i, rhomt[ia][i])
    end

    #exit()

    # interstitial density determined from the atomic tails and excess charge
    println("zfft before = ", sum(zfft))
    G_to_R!(pw, zfft)
    zfft[:] = zfft[:]*Npoints
    println("zfft after = ", sum(zfft))
    for ip in 1:Npoints
        rhoir[ip] = real(zfft[ip]) + t1
        # make sure that the density is always positive
        if rhoir[ip] < 1.e-10
            rhoir[ip] = 1.e-10
        end
    end

    println("Npoints = ", Npoints)

    @printf("sum rhoir = %18.10e\n", sum(rhoir))
    ss = 0.0
    for ia in 1:Natoms
        ss = ss + sum(rhomt[ia])
    end
    @printf("sum rhomt = %18.10e\n", ss)
    @printf("Total = %18.10e\n", ss + sum(rhoir))

    return
end