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

    nrsp = atsp_vars.nrsp
    rsp = atsp_vars.rsp
    rhosp = atsp_vars.rhosp 

    lmaxi = mt_vars.lmaxi
    nrmt = mt_vars.nrmt
    npmt = mt_vars.npmt

    sfacg = calc_strfact(atoms, pw)

    lmax = min(lmaxi,1)

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
    ffg = zeros(Float64,Ng,Nspecies)

    # ALLOCATE(ffg(ngvec), wr(nrspmax), fr(nrspmax))
    for isp in 1:Nspecies

        # local arrays inside the loop
        wr = zeros(Float64,nrsp[isp])
        fr = zeros(Float64,nrsp[isp])

        nr = nrmt[isp]
        nrs = nrsp[isp]
        nro = nrs - nr + 1

        # determine the weights for the radial integral
        @views wsplint!( nro, rsp[isp][nr:nr+nro-1] ,wr[nr:nr+nro-1] )
        for ig in 1:Ng
            t1 = sqrt(G2[ig])
            # spherical bessel function j_0(x) times the atomic density tail
            if t1 > epslat
                t2 = 1.0/t1
                for ir in nr:nrs
                  x = t1*rsp[isp][ir]
                  fr[ir] = t2*sin(x)*rhosp[isp][ir]*rsp[isp][ir]
                end
            else
                for ir in nr:nrs
                    fr[ir] = rhosp[isp][ir] * rsp[isp][ir]^2
                end
            end
            @views t1 = dot( wr[nr:nrs], fr[nr:nrs] )
            # apply low-pass filter
            t1 = t1*exp(-4.0*G2[ig]/gmaxvr^2)
            ffg[ig,isp] = (4*pi/CellVolume)*t1
        end
    end
        
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ig in 1:Ng
            ip = pw.gvec.idx_g2r[ig]
            #zfft[ip] = zfft[ip] + ffg[ig,isp]*conj(sfacg[ig,isp])
            zfft[ip] = zfft[ip] + ffg[ig,isp]*sfacg[ig,isp]
        end
    end

    println("sum zfft = ", sum(zfft))
    #DEALLOCATE(ffg, wr, fr)


    return
end