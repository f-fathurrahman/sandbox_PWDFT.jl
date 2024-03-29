function electrons_scf_broyden!(
    Ham::Hamiltonian, psiks;
    NiterMax=150,
    betamix=0.7,
    etot_conv_thr=1e-6,
    ethr_evals_last=1e-13
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

    evals = zeros(Float64, Nstates, Nkpt*Nspin)
    
    mixdim = 8
    
    # Mix directly in R-space
    df = zeros(Float64,Npoints*Nspin, mixdim)
    dv = zeros(Float64,Npoints*Nspin, mixdim)
    
    #df = zeros(ComplexF64,Npoints*Nspin, mixdim)
    #dv = zeros(ComplexF64,Npoints*Nspin, mixdim)

    #Ng = Ham.pw.gvec.Ng
    #df = zeros(ComplexF64, Ng*Nspin, mixdim)
    #dv = zeros(ComplexF64, Ng*Nspin, mixdim)

    ethr = 1e-5 # default

    is_converged = false

    for iterSCF in 1:NiterMax
        
        # determine convergence criteria for diagonalization
        if iterSCF == 1
            ethr = 0.1
        elseif iterSCF == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ethr_evals_last )
        end

        println("\niterSCF = ", iterSCF)
        println("Davidson diagonalization with ethr = ", ethr)
        evals[:,:] = diag_davidson_qe!( Ham, psiks, tol=ethr )

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
        #mix_broyden_03!( Ham, Rhoe, Rhoe_new, betamix, iterSCF, mixdim, df, dv )
        #mix_broyden_02!( Ham.pw, Rhoe, Rhoe_new, betamix, iterSCF, mixdim, df, dv )
        
        mix_broyden_02!( Rhoe, Rhoe_new, betamix, iterSCF, mixdim, df, dv )
        #mix_anderson!( Rhoe, Rhoe_new, betamix, df, dv, iterSCF, mixdim )

        println("integ Rhoe after mix: ", sum(Rhoe)*dVol)

        diffRhoe = norm(Rhoe - Rhoe_new)
        #diffRhoe = dot(Rhoe - Rhoe_new, Rhoe - Rhoe_new)

        #
        Ehartree, Exc, Evtxc = update_from_rhoe!(Ham, Rhoe)

        descf = sum( (Rhoe[:,1] .- Rhoe_new[:,1]).*(Vhartree + Vxc[:,1]) )*dVol

        Etot = Eband + deband + Ehartree + Exc + Ham.energies.NN + descf # entropy missed
        
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
        @printf("Total    = %18.10f\n", 2*Etot)

        #
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

        #if diffRhoe < 1e-6
        #    @printf("converged by diffRhoe\n")
        #    is_converged = true
        #    return
        #end

        Etot_old = Etot
        flush(stdout)
    end

    if !is_converged
        @printf("WARNING: SCF is not converged after %d iterations\n", NiterMax)
    end

    # Evaluating energy
    energies = calc_energies(Ham, psiks)
    println("\nUsing original formula for total energy")
    println(energies)

    return
end