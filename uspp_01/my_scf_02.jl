function my_scf!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.2,
    etot_conv_thr=1e-6
)

    Ham.energies.NN = calc_E_NN(Ham.atoms)

    Etot_old = 0.0
    Nconv = 0

    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    dVol = CellVolume/Npoints
    Nstates = Ham.electrons.Nstates

    Rhoe = Ham.rhoe
    println("Initial integ Rhoe = ", sum(Rhoe)*dVol)

    diffRhoe = 0.0

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")

    Vin = zeros(Float64, Npoints)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    Vhartree = Ham.potentials.Hartree
    Vxc = Ham.potentials.XC
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk

    mixer = BroydenMixer(Rhoe, betamix, mixdim=8)

    evals = Ham.electrons.ebands
    is_converged = false
    
    # Other energy terms
    Eband = 0.0
    deband = 0.0
    descf = 0.0
    Etot = 0.0
    Ehartree = 0.0
    Exc = 0.0

    for iterSCF in 1:NiterMax
        
        println("\niterSCF = ", iterSCF)
        evals[:,:] .= diag_davidson_qe!( Ham, psiks )

        println("Eigenvalues in eV: ikspin = 1")
        for ist in 1:Nstates
            @printf("%5d %18.10f\n", ist, evals[ist,1]*Ha2eV)
        end

        Eband = 0.0
        for ispin in 1:Nspin, ik in 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            for ist in 1:Nstates
                Eband = Eband + wk[ik]*Focc[ist,ikspin]*evals[ist,ikspin]
            end
        end

        Rhoe_new = calc_rhoe_uspp( Ham, psiks )
        println("integ Rhoe_new = ", sum(Rhoe_new)*dVol)

        # 
        Vin[:] = Vhartree[:] + Vxc[:,1]
        deband = -sum(Vin .* Rhoe_new[:,1])*dVol

        #
        # Mix the density
        #
        #Rhoe = betamix*Rhoe_new + (1 - betamix)*Rhoe
        do_mix!(mixer, Rhoe, Rhoe_new, iterSCF)

        println("integ Rhoe after mix: ", sum(Rhoe)*dVol)

        diffRhoe = norm(Rhoe - Rhoe_new)

        #
        Ehartree, Exc, Evtxc = update_from_rhoe!(Ham, Rhoe)

        descf = sum( (Rhoe[:,1] .- Rhoe_new[:,1]).*(Vhartree + Vxc[:,1]) )*dVol

        # FIXME: Entropy term is missing (metallic occupation is not yet supported)
        Etot = Eband + deband + Ehartree + Exc + Ham.energies.NN + descf # entropy missed

        diffEtot = abs(Etot - Etot_old)
        @printf("\nSCF: %5d %18.10f %10.5e %10.5e\n", iterSCF, Etot, diffEtot, diffRhoe)
        if diffEtot <= etot_conv_thr
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        if Nconv >= 2
            @printf("SCF is converged in %d iterations\n", iterSCF)
            is_converged = true
            break
        end

        Etot_old = Etot
        flush(stdout)
    end

    if !is_converged
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    println()
    println("Final result")
    println()
    println("-----------------------")        
    println("Energy components in Ry")
    println("-----------------------")
    @printf("Eband    = %18.10f\n", Eband*2)
    @printf("deband   = %18.10f\n", deband*2)
    @printf("descf    = %18.10f\n", descf*2)
    @printf("-----------------------------\n")
    @printf("OneEle   = %18.10f\n", 2*(Eband + deband))
    @printf("Ehartree = %18.10f\n", 2*Ehartree)
    @printf("Exc      = %18.10f\n", 2*Exc)
    @printf("NN       = %18.10f\n", 2*Ham.energies.NN)
    @printf("-----------------------------\n")
    @printf("! Total  = %18.10f\n", 2*Etot)

    return
end
