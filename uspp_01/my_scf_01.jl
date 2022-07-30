function my_scf!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.1,
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

    for iterSCF in 1:NiterMax
        
        println("\niterSCF = ", iterSCF)
        evals = diag_davidson_qe!( Ham, psiks )

        println("Eigenvalues: ")
        for ist in 1:Nstates
            @printf("%5d %18.10f\n", ist, evals[ist])
        end

        println("Ortho check after diag_davidson_qe")
        ortho_check_with_S(Ham, psiks[1])

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
            @printf("SCF is converged in %d iterations\n", iterSCF)
            return
        end

        if diffRhoe < 1e-5
            @printf("converged by diffRhoe\n")
            return
        end

        Etot_old = Etot
        flush(stdout)
    end
    @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    return
end