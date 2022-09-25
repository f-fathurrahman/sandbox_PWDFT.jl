function my_scf!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.2,
    etot_conv_thr=1e-6
)

    Ham.energies.NN = calc_E_NN(Ham.atoms)

    Etot_old = 0.0
    Nconv = 0

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)
    Nstates = Ham.electrons.Nstates

    Rhoe = Ham.rhoe
    println("Initial integ Rhoe = ", sum(Rhoe)*dVol)

    diffRhoe = 0.0

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")

    is_converged = false

    evals = Ham.electrons.ebands

    for iterSCF in 1:NiterMax
        
        println("\niterSCF = ", iterSCF)
        evals[:,:] .= diag_davidson_qe!( Ham, psiks )

        ikspin = 1
        println("Eigenvalues in eV: ")
        for ist in 1:Nstates
            @printf("%5d %18.10f\n", ist, evals[ist,ikspin]*Ha2eV)
        end

        Rhoe_new = calc_rhoe_uspp( Ham, psiks )
        println("integ Rhoe_new = ", sum(Rhoe_new)*dVol)

        Rhoe = betamix*Rhoe_new + (1 - betamix)*Rhoe
        println("integ Rhoe after mix: ", sum(Rhoe)*dVol)

        diffRhoe = norm(Rhoe - Rhoe_new)

        #
        _ = update_from_rhoe!(Ham, Rhoe)
        #
        Ham.energies = calc_energies(Ham, psiks)
        println(Ham.energies)
        #
        Etot = sum(Ham.energies)
        #
        diffEtot = abs(Etot - Etot_old)
        @printf("SCF: %5d %18.10f %10.5e %10.5e\n", iterSCF, Etot, diffEtot, diffRhoe)
        if diffEtot <= etot_conv_thr
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        if Nconv >= 2
            is_converged = true
            @printf("SCF is converged in %d iterations\n", iterSCF)
            break
        end

        Etot_old = Etot
        flush(stdout)
    end

    if !is_converged
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    energies = Ham.energies

    E_OneEle = energies.Kinetic + energies.Ps_loc + energies.Ps_nloc

    println("\nTotal energy in Ry")
    println("-------------------------------------")
    @printf("Kinetic    energy: %18.10f\n", 2*energies.Kinetic )
    @printf("Ps_loc     energy: %18.10f\n", 2*energies.Ps_loc )
    @printf("Ps_nloc    energy: %18.10f\n", 2*energies.Ps_nloc )
    @printf("OneEle     energy: %18.10f\n", 2*E_OneEle) 
    @printf("Hartree    energy: %18.10f\n", 2*energies.Hartree )
    @printf("XC         energy: %18.10f\n", 2*energies.XC )
    @printf("-------------------------------------\n")
    
    E_elec = energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
             energies.Hartree + energies.XC + energies.mTS
    
    @printf("Electronic energy: %18.10f\n", 2*E_elec)
    @printf("NN         energy: %18.10f\n", 2*energies.NN )
    @printf("-------------------------------------\n")
    
    E_total = E_elec + energies.NN
    @printf("Total      energy: %18.10f\n", 2*E_total)

    return
end
