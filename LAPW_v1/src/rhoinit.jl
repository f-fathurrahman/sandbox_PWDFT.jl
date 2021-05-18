function rhoinit!(
    atoms, atsp_vars,
    mt_vars,
    pw
)

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

    sfacg = calc_strfact(atoms, pw)
    ylmg = zeros(ComplexF64, lmmaxo, Ng)
    genylmg!(mt_vars.lmaxo, pw.gvec.G, ylmg)

    lmax = min(mt_vars.lmaxi,1)

    epslat = 1e-6
  
    # zero the charge density arrays
    rhomt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        rhomt[ia] = zeros(Float64, npmt[isp])
    end
    rhoir = zeros(Float64,Npoints)
  
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

        println("nr  = ", nr)
        println("nrs = ", nrs)
        println("nro = ", nro)

        # determine the weights for the radial integral
        # integrate for radial points outside the muffin-tin
        idx = nr:nr+nro-1 # or better: idx = nr:nrs
        println("idx = ", idx)
        @views wsplint!( nro, rsp[isp][idx] ,wr[idx] )
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

    println("sum sfacg = ", sum(sfacg))
    s = 0.0
    for isp in 1:Nspecies
        s = s + sum(rhosp[isp])
    end
    println("sum rhosp = ", s)
    println("sum zfft = ", sum(zfft))

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
        println("nrci = ", nrci)
        println("nrc  = ", nrc)
        zfmt[1:npcmt[isp]] .= 0.0
        for ig in 1:Ng
            ip = idx_g2r[ig]
            for irc in 1:nrc
                x = sqrt(G2[ig])*rcmt[isp][irc]
                @views sbessel!(lmax, x, jl[:,irc])
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
        z_to_rf_mt!( mt_vars, nrc, nrci, zfmt, rhomt[isp] )
    end

    # convert the density from a coarse to a fine radial mesh
    rf_mt_c_to_f!( atoms, atsp_vars, mt_vars, rhomt )

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
            rhomt[isp][i] = rhomt[isp][i] + t2
            i = i + lmmaxi
        end
        for ir in nri+1:nr
            t2 = (t1 + rhosp[isp][ir])/y00
            rhomt[isp][i] = rhomt[isp][i] + t2
            i = i + lmmaxo
        end
    end

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
    for isp in 1:Nspecies
        ss = ss + sum(rhomt[isp])
    end
    @printf("sum rhomt = %18.10e\n", ss)
    @printf("Total = %18.10e\n", ss + sum(rhoir))

    return
end