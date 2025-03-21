function genlofr!(atoms, eqatoms, mt_vars, apwlo_vars, vsmt)

    y00 = 0.5/sqrt(pi)

    Natoms = atoms.Natoms
    done = zeros(Bool, Natoms)
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    nrmt = mt_vars.nrmt
    nrmti = mt_vars.nrmti
    nrmtmax = maximum(nrmt)
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo
    rlmt = mt_vars.rlmt
    rmt = mt_vars.rmt

    lorbordmax = apwlo_vars.lorbordmax
    deapwlo = apwlo_vars.deapwlo
    lorbord = apwlo_vars.lorbord
    lorbe = apwlo_vars.lorbe
    lorbdm = apwlo_vars.lorbdm
    nlorb = apwlo_vars.nlorb
    lorbl = apwlo_vars.lorbl
    nplorb = apwlo_vars.nplorb
    nlomax = apwlo_vars.nlomax
    lofr = apwlo_vars.lofr

    p0 = zeros(Float64, nrmtmax, lorbordmax)
    p1 = zeros(Float64, nrmtmax)
    q0 = zeros(Float64, nrmtmax)
    q1 = zeros(Float64, nrmtmax)
    ep0 = zeros(Float64, nrmtmax, lorbordmax)
    p0s = zeros(Float64, nrmtmax, nlomax)
    ep0s = zeros(Float64, nrmtmax, nlomax)
    vr = zeros(Float64, nrmtmax)
    fr = zeros(Float64, nrmtmax)

    xa = zeros(Float64, nplorb)
    ya = zeros(Float64, nplorb)
    A = zeros(Float64, nplorb, nplorb)
    b = zeros(Float64, nplorb)

    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = nrmt[isp]
        nri = nrmti[isp]
        if done[ia]
            continue
        end
        # use spherical part of potential
        i = 1
        for ir in 1:nri
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxi
        end
        for ir in (nri+1):nr
            vr[ir] = vsmt[ia][i]*y00
            i = i + lmmaxo
        end
        @views rgrid = rlmt[isp][:,1] # nr will be taken from this
        #
        for ilo in 1:nlorb[isp]
            l = lorbl[isp][ilo]
            for jo in 1:lorbord[isp][ilo]
                # linearisation energy accounting for energy derivative
                E = lorbe[ia][ilo][jo] + lorbdm[isp][ilo][jo]*deapwlo
                # integrate the radial Schrodinger equation
                @views p0view = p0[1:nr,jo]
                _ = rschrodint!(l, E, rgrid, vr, p0view, p1, q0, q1)
                #@printf("\nl, ilo, linearization energy = %4d %4d %18.10f %4d\n", l, ilo, E, nn)
                ep0[1:nr,jo] .= E*p0[1:nr,jo]
                # normalise radial functions
                fr[1:nr] = p0[1:nr,jo].^2
                t1 = splint(nr, rgrid, fr)
                t1 = 1.0/sqrt(abs(t1))
                p0[1:nr,jo] .*= t1
                ep0[1:nr,jo] .*= t1
                #@printf("jo, sum abs ep0 = %4d %18.10f\n", jo, sum(abs.(ep0[1:nr,jo])))
                # set up the matrix of radial derivatives
                for i in 1:nplorb
                    ir = nr - nplorb + i
                    xa[i] = rlmt[isp][ir,1]
                    ya[i] = p0[ir,jo]*rlmt[isp][ir,-1]
                end
                for io in 1:lorbord[isp][ilo]
                    A[io,jo] = polynm(io-1, nplorb, xa, ya, rmt[isp])
                end 
            end 
            # set up the target vector
            b[:] .= 0.0 # zero
            b[lorbord[isp][ilo]] = 1.0 # except this one
            #CALL dgesv(lorbord(ilo,is), 1, a, nplorb, ipiv, b, nplorb, info)
            idx1 = 1:lorbord[isp][ilo]
            b[idx1] = A[idx1,idx1]\b[idx1]
            #IF(info != 0) goto 10
            # generate linear superposition of radial functions
            p0s[1:nr,ilo] .= 0.0
            ep0s[1:nr,ilo] .= 0.0
            for io in 1:lorbord[isp][ilo]
                t1 = b[io]
                #@printf("io = %4d b = %18.10f\n", io, b[io])
                p0s[1:nr,ilo] .+= t1*p0[1:nr,io]
                ep0s[1:nr,ilo] .+= t1*ep0[1:nr,io]
            end 
            # normalize radial functions
            fr[1:nr] .= p0s[1:nr,ilo].^2
            t1 = splint(nr, rgrid, fr)
            t1 = 1.0/sqrt(abs(t1))
            #@printf("Scaling factor t1 from p0s = %18.10f\n", t1)
            p0s[1:nr,ilo] .*= t1
            ep0s[1:nr,ilo] .*= t1
            #@printf("ilo, sum abs ep0s = %4d %18.10f\n", ilo, sum(abs.(ep0s[1:nr,ilo])))
            #
            # subtract linear combination of previous local-orbitals with same l
            for jlo in 1:(ilo-1)
                if lorbl[isp][jlo] == l
                    fr[1:nr] .= p0s[1:nr,ilo] .* p0s[1:nr,jlo]
                    t1 = -splint(nr, rgrid, fr)
                    p0s[1:nr,ilo] .+= t1*p0s[1:nr,jlo]
                    ep0s[1:nr,ilo] .+= t1*ep0s[1:nr,jlo]
                end
            end
            # normalize radial functions again
            fr[1:nr] = p0s[1:nr,ilo].^2
            t1 = splint(nr, rgrid, fr)
            t1 = abs(t1)
            if t1 < 1e-25
                error("Degenerate LO radial functions")
            end
            t1 = 1.0/sqrt(t1)
            #@printf("Scaling factor t1 again = %18.10f\n", t1)
            p0s[1:nr,ilo] .*= t1
            ep0s[1:nr,ilo] .*= t1
            # divide by r and store in global array
            for ir in 1:nr
                t1 = rlmt[isp][ir,-1]
                lofr[ia][ilo][ir,1] = t1*p0s[ir,ilo]
                lofr[ia][ilo][ir,2] = t1*ep0s[ir,ilo]
            end
        end # for ilo
        done[ia] = true
        # copy to equivalent atoms
        for ja in 1:Natoms
            if !done[ja] && eqatoms[ia,ja]
                for ilo in 1:nlorb[isp]
                    lofr[ja][ilo][1:nr,1] .= lofr[ia][ilo][1:nr,1]
                    lofr[ja][ilo][1:nr,2] .= lofr[ia][ilo][1:nr,2]
                end 
                done[ja] = true
            end
        end
    end
    return
end
