# Need to run setup first


function main_cg_01(Ham; NiterMax=100, psiks=nothing, Haux=nothing)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Rhoe = Ham.rhoe
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    # Initialize electronic variables: `psiks` and `Haux`:
    Random.seed!(1234)
    if isnothing(psiks)
        psiks = rand_BlochWavefunc(Ham)
    end

    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end
    # Calculate Hsub
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        Hsub[ikspin][:,:] = psiks[ikspin]' * (Ham * psiks[ikspin])
    end

    # Prepare Haux (random numbers)
    #
    if isnothing(Haux)
        Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
        for ikspin in 1:Nkspin
            Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
            # the same as Hsub
            Haux[ikspin][:,:] = Hsub[ikspin][:,:]
            Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
        end
    end

    # Gradients, subspace Hamiltonian
    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    d = zeros_BlochWavefunc(Ham)
    gPrev = zeros_BlochWavefunc(Ham)
    #
    g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    d_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    gPrev_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        d_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        gPrev_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end

    rots_cache = RotationsCache(Nkspin, Nstates);

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache )

    # Update Hamiltonian, compute energy and gradients at current psiks and Haux:

    # Update Hamiltonian before evaluating free energy
    update_from_ebands!( Ham )
    #update_from_wavefunc!( Ham, psiks )
    E1 = calc_Lfunc( Ham, psiks )
    println("E1 = $(E1)")
    #
    # Calculate gradients
    calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    println("Initial Focc = ")
    display(Ham.electrons.Focc); println()
    println("Initial ebands (w.r.t Fermi) = ")
    display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()

    α_t_start = 1.0
    α_t_min = 1e-5
    α_t = α_t_start

    do_force_grad_dir = true
    gKNorm = 0.0
    gKNormPrev = 0.0
    # current and previous norms of the preconditioned gradient

    ok_paw = any(Ham.pspotNL.are_paw) # only used for printing?

    for iterCG in 1:NiterMax

        println("\nStart iterCG = ", iterCG)

        gKNorm = 2*real(dot(g, Kg)) + real(dot(g_Haux, Kg_Haux))

        β = 0.0
        if !do_force_grad_dir
            gPrevKg = 2*real(dot(gPrev, Kg)) + real(dot(gPrev_Haux, Kg_Haux))
            gd = 2*real(dot(g, d)) + real(dot(g_Haux, d_Haux))
            gg = 2*real(dot(g, g)) + real(dot(g_Haux, g_Haux))
            dd = 2*real(dot(d, d)) + real(dot(d_Haux, d_Haux))
            if gg*dd > 0
                @printf("linmin: %10.3le\n", gd/sqrt(gg*dd))
            else
                @warn "Negative gg*dd encountered"
            end
            if gKNorm*gKNormPrev > 0
                @printf("cgtest: %10.3le\n", gPrevKg/sqrt(gKNorm*gKNormPrev))
            else
                @warn "Negative gKNorm*gKNormPrev encountered"
            end
            # Update beta:
            println("gKNorm = $(gKNorm), gPrevKg = $(gPrevKg)")
            β = (gKNorm - gPrevKg)/gKNormPrev
            println("β = ", β)              
            if β < 0.0
                println("!!!! Resetting CG because β is negative")
                β = 0.0
            end
        end

        do_force_grad_dir = false

        # XXX TODO Check convergence here?

        # Save previous gradient
        gKNormPrev = gKNorm
        for ikspin in 1:Nkspin
            gPrev[ikspin][:] = g[ikspin][:]
            gPrev_Haux[ikspin][:] = g_Haux[ikspin][:]
        end

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin] = -Kg[ikspin] + β*d[ikspin]
            d_Haux[ikspin] = -Kg_Haux[ikspin] + β*d_Haux[ikspin]
        end
        constrain_search_dir!(Ham, d, psiks)

        #
        # Do line minimization:
        E_new, is_success, α = linmin_quad_v01!(
            α_t,
            Ham, psiks, Haux, Hsub, g, g_Haux, Kg, Kg_Haux, d, d_Haux, rots_cache, E1
        )
        println("Test grad psiks before rotate: $(2*dot(g, psiks))")
        println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
        rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
 
        #
        if is_success
            α_t = α
            println("linminQuad is successful. α_t is updated to α = ", α)
            if α_t < α_t_min
                # bad step size: make sure next test step size is not too bad
                α_t = α_t_start 
                println("Bad step size is encountered, α_t is set to α_t_start = \n", α_t_start)
            end
        else
            println("WARN: Line minimization is not successful")
            #
            do_step_psiks_Haux!(-α, Ham, psiks, Haux, d, d_Haux, rots_cache)
            # calculate energy and gradients
            update_from_ebands!( Ham )
            update_from_wavefunc!( Ham, psiks )
            # Now, we are ready to evaluate
            E_new = Inf # calc_Lfunc( Ham, psiks )
            calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
            calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
            #
            if β > 0.0
                # Failed, but not along the gradient direction:
                println("Forcing gradient direction")
                do_force_grad_dir = true
            else
                # Failed along the gradient direction
                println("Probably round off error")
                break
            end
        end

        if Nspin == 2
            magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe integ magn = $magn")
        else
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe")
        end
        ΔE = abs(E_new - E1)
        println("iterCG: $(iterCG) E_new = $(E_new) ΔE = $(ΔE)")
        println("Focc = ")
        display(Ham.electrons.Focc); println()
        println("ebands (w.r.t) Fermi energy = ")
        display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()
        println("Energies:")
        println(Ham.energies, use_smearing=true, is_paw=ok_paw)

        if ΔE < 1e-8
            println("\nConverged !!!")
            break
        end 
        
        # New iterations, variables are updated in linmin_quad_v01
        E1 = E_new
    end

    serialize("psiks.jldat", psiks)
    serialize("Haux.jldat", Haux)

    #@infiltrate

    return
end

