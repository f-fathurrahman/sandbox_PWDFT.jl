function my_calc_forces_scf_corr!(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    potentials::Potentials,
    F_scf_corr::Matrix{Float64}
)

    # Calculate difference in potential
    Vtot = potentials.Total
    VtotOld = potentials.TotalOld
    
    eps8 = 1e-8
    G2_shells = pw.gvec.G2_shells
    idx_g2r = pw.gvec.idx_g2r
    G = pw.gvec.G
    idx_g2shells = pw.gvec.idx_g2shells
    Ngl = length(G2_shells)
    Ng = pw.gvec.Ng

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    X = atoms.positions

    Nspin = size(Vtot,2)
    Npoints = size(Vtot,1)

    println("Npoints = ", Npoints)
    println("pw.Ns = ", pw.Ns)

    # Calculate difference between new and old potential
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin, ip in 1:Npoints
        ctmp[ip] += Vtot[ip,ispin] - VtotOld[ip,ispin]
    end

    println("sum ctmp = ", sum(ctmp))

    #R_to_G!(pw, ctmp)
    ff = reshape(ctmp, pw.Ns)
    # Use FFT again
    planfw = plan_fft!( zeros(ComplexF64,pw.Ns) ) # using default plan
    planfw*ff

    ctmp[:] *= (1/Npoints) # rescale

    println("pw.planfw = ", typeof(pw.planfw))

    println("sum ctmp after forward FFT = ", sum(ctmp))

    # Determine maximum radial points for every atomic species
    NpointsMax = 0
    for isp in 1:Nspecies
        if NpointsMax < pspots[isp].Nr
            NpointsMax = pspots[isp].Nr
        end
    end
    aux = zeros(Float64, NpointsMax)
    rhocgnt = zeros(Float64, Ngl)

    fill!(F_scf_corr, 0.0)

    for isp in 1:Nspecies
        psp = pspots[isp]
        # G != 0 terms
        for igl in 2:Ngl
            gx = sqrt(G2_shells[igl])
            for ir in 1:psp.Nr
                if psp.r[ir] < eps8
                   aux[ir] = psp.rhoatom[ir]
                else
                   aux[ir] = psp.rhoatom[ir]*sin(gx*psp.r[ir])/(psp.r[ir]*gx)
                end
            end
            rhocgnt[igl] = PWDFT.integ_simpson( psp.Nr, aux, psp.rab )
        end

        # sum over atoms
        for ia in 1:Natoms
            if isp != atm2species[isp]
                continue
            end
            #
            for ig in 2:Ng
                igl = idx_g2shells[ig]
                ip = idx_g2r[ig]
                GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[3,ig]*X[3,ia]
                Sf = sin(GX) + im*cos(GX)
                @views F_scf_corr[:,ia] .+= real.( rhocgnt[igl] * Sf * G[:,ig] * conj(ctmp[ip]) )
            end
        end

    end

    return

end