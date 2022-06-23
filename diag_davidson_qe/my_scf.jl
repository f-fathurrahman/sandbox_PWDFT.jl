function my_scf!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    NiterMax=150,
    betamix=0.2,
    etot_conv_thr=1e-6
)
    
    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, Rhoe)
    
    Rhoe_new = similar(Rhoe)
    
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    Etot_old = 0.0
    Nconv = 0

    @printf("\n")
    @printf("SCF iteration starts (with density mixing), betamix = %f\n", betamix)
    @printf("\n")

    for iterSCF in 1:NiterMax
        #_ = diag_davidson_qe!( Ham, psiks )
        _ = diag_LOBPCG!( Ham, psiks )
        Rhoe_new = calc_rhoe( Ham, psiks )
        Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe
        update!(Ham, Rhoe)
        Ham.energies = calc_energies(Ham, psiks)
        Etot = sum(Ham.energies)
        diffEtot = abs(Etot - Etot_old)
        @printf("%5d %18.10f %18.10e\n", iterSCF, Etot, diffEtot)
        if diffEtot <= etot_conv_thr
            Nconv = Nconv + 1
        else
            Nconv = 0
        end
        if Nconv >= 2
            @printf("SCF is converged in %d iterations\n", iterSCF)
            return
        end
        Etot_old = Etot
        flush(stdout)
    end
    @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    return
end