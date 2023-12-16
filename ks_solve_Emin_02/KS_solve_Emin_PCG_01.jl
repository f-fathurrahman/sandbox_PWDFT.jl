# also include overlap operator
function new_calc_grad( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    
    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt
    #
    Focc = Ham.electrons.Focc
    #
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]
    #
    grad = zeros( ComplexF64, Ngw, Nstates )

    H_psi = op_H( Ham, psi )
    S_psi = op_S( Ham, psi )
    for ist in 1:Nstates
        grad[:,ist] = H_psi[:,ist]
        for jst in 1:Nstates
            grad[:,ist] = grad[:,ist] - dot( psi[:,jst], H_psi[:,ist] ) * S_psi[:,jst]
        end
        grad[:,ist] = Focc[ist,ikspin]*grad[:,ist]*Ham.pw.gvecw.kpoints.wk[ik]
    end

    # the usual case of constant occupation numbers
    if all(Focc[:,ikspin] .== 2.0) == true || all(Focc[:,ikspin] .== 1.0) == true
        # immediate return
        return grad
    end

#    F = Matrix(Diagonal(Focc[:,ikspin]))
#    ℍ = psi' * H_psi
#    HFH = ℍ*F - F*ℍ
#    ℚ = 0.5*HFH
#    grad[:,:] = grad[:,:] + psi*ℚ

    return grad

end





function KS_solve_Emin_PCG_01!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    startingrhoe::Symbol=:gaussian,
    skip_initial_diag::Bool=false,
    α_t::Float64=3e-5,
    NiterMax::Int64=200,
    verbose::Bool=true,
    print_final_ebands::Bool=false,
    print_final_energies::Bool=true,
    i_cg_beta::Int64=2,
    etot_conv_thr::Float64=1e-6,
    α_max::Float64=2.1,
    restrict_linmin::Bool=false
)

    Rhoe = Ham.rhoe # alias
    PWDFT._prepare_scf!(Ham, psiks)
    # Can be used in Emin too, despite the name

    pw = Ham.pw
    electrons = Ham.electrons
    
    Focc = electrons.Focc
    Nstates = electrons.Nstates
    Nelectrons = electrons.Nelectrons
    
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Nkpt = pw.gvecw.kpoints.Nkpt
    
    Nspin = electrons.Nspin
    Nkspin = Nkpt*Nspin

    #
    # Variables for PCG
    #
    g      = zeros_BlochWavefunc( Ham )
    d      = deepcopy(g)
    g_old  = deepcopy(g)
    d_old  = deepcopy(g)
    Kg     = deepcopy(g)
    Kg_old = deepcopy(g)
    psic   = deepcopy(g)
    gt     = deepcopy(g)
    
    β = zeros(Nkspin)
    α = zeros(Nkspin)

    Etot_old = 0.0

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # Calculate energy at this psi
    energies = calc_energies( Ham, psiks )
    Ham.energies = energies
    Etot = sum(energies)

    ok_paw = any(Ham.pspotNL.are_paw)

    Nconverges = 0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG\n")
        @printf("-------------------------------------\n")
        @printf("NiterMax  = %d\n", NiterMax)
        @printf("α_t       = %e\n", α_t)
        @printf("conv_thr  = %e\n", etot_conv_thr)
        @printf("\n")
    end


    for iter in 1:NiterMax

        for ispin in 1:Nspin, ik in 1:Nkpt

            Ham.ik = ik
            Ham.ispin = ispin
            i = ik + (ispin - 1)*Nkpt

            g[i] = new_calc_grad( Ham, psiks[i] )
            Kg[i] = Kprec( Ham.ik, pw, g[i] )

            # XXX: define function trace for real(sum(conj(...)))
            if iter != 1
                β[i] = real( dot(g[i]-g_old[i], Kg[i]) )/real( dot(g_old[i],Kg_old[i]) )
            end
            if β[i] < 0.0 || isnan(β[i])
                #println("Resetting β")
                β[i] = 0.0
            end

            d[i] = -Kg[i] + β[i] * d_old[i]

            psic[i] = psiks[i] + α_t*d[i]
            ortho_sqrt!(Ham, psic[i])
        end # ik, ispin
        
        calc_rhoe!( Ham, psic, Rhoe )

        update_from_rhoe!(Ham, psic, Rhoe)

        for ispin in 1:Nspin, ik in 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            i = ik + (ispin - 1)*Nkpt
            gt[i] = new_calc_grad(Ham, psic[i])
            denum = real(sum(conj(g[i]-gt[i]).*d[i]))
            if denum != 0.0
                α[i] = abs( α_t*real(sum(conj(g[i]).*d[i]))/denum )
            else
                α[i] = 0.0
            end

            if restrict_linmin
                if α[i] > α_max
                    @printf("α for ikspin #%d is too large: %f, restrict it to %f\n", i, α[i], α_max)
                    α[i] = α_max
                end
            end

            # Update wavefunction
            psiks[i] = psiks[i] + α[i]*d[i]
            ortho_sqrt!(Ham, psiks[i])
        end

        # Update Rhoe and potentials
        calc_rhoe!( Ham, psiks, Rhoe )
        update_from_rhoe!(Ham, psiks, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diffE = Etot_old - Etot

        if verbose
            @printf("Emin_PCG step %4d = %18.10f %16.6e\n", iter, Etot, diffE)
        end
        
        if verbose && (diffE < 0.0) && (iter > 1)
            println("*** WARNING: Etot is not decreasing")
        end
        
        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            if verbose
                @printf("\nEmin_PCG is converged in iter: %d\n", iter)
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot

        flush(stdout)
    end

    if verbose && print_final_energies
        @printf("\n")
        @printf("-------------------------\n")
        @printf("Final Kohn-Sham energies:\n")
        @printf("-------------------------\n")
        @printf("\n")
        println(Ham.energies, is_paw=ok_paw)
    end

    return

end
