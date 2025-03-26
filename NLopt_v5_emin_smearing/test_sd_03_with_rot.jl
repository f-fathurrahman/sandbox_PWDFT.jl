using Revise
using LinearAlgebra
using Printf
using Infiltrator
using Random
using PWDFT

includet("smearing.jl")
includet("occupations.jl")
includet("Lfunc.jl")
includet("gradients_psiks_Haux.jl")
includet("utilities_emin_smearing.jl")
includet("prepare_Ham_various.jl")


function main_sd_03(Ham; NiterMax=100)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates;

    # Initialize electronic variables: `psiks` and `Haux`:
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham);

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
    #@infiltrate

    # Prepare Haux (random numbers)
    #
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    # For Haux, choose between generic symmetric Haux:
    #=
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
        # the same as Hsub
        Haux[ikspin][:,:] = Hsub[ikspin][:,:]
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    end
    =#

    # eigenvalues of Hsub
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm(0 => eigvals(Hsub[ikspin]))
    end

    # or diagonal Haux:
    # Prepare Haux (random numbers)
    #=
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm(0 => sort(randn(Float64, Nstates)))
    end
    =#

    # Gradients, subspace Hamiltonian
    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    d = zeros_BlochWavefunc(Ham)
    #
    g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    d_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        d_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end

    rots_cache = RotationsCache(Nkspin, Nstates);
    Haux_orig = copy(Haux);

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache,
        do_ortho_psi=false, overwrite_Haux=true
    )

    # Update Hamiltonian, compute energy and gradients at current psiks and Haux:

    # Update Hamiltonian before evaluating free energy
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks )
    E1 = calc_Lfunc( Ham, psiks )
    println("E1 = $(E1)")
    #
    # Calculate gradients
    calc_grad_psiks!(Ham, psiks, g, Hsub)
    my_Kprec!(Ham, g, Kg)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    println("Initial Focc = ")
    display(Ham.electrons.Focc); println
    println("Initial ebands (w.r.t Fermi) = ")
    display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println

    α_t_start = 1.0
    α_t_min = 1e-10
    α_t = α_t_start

    Rhoe = Ham.rhoe
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    for iterSD in 1:NiterMax

        println("\nStart iterSD = ", iterSD)

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin][:,:] = -Kg[ikspin][:,:]
            d_Haux[ikspin][:,:] = -Kg_Haux[ikspin][:,:]
        end
        constrain_search_dir!(d, psiks)

        reset_rotations!(rots_cache)

        #
        # Do line minimization:
        E_new, is_success, α = linmin_quad_v01!(
            α_t,
            Ham, psiks, Haux, Hsub, g, g_Haux, Kg, Kg_Haux, d, d_Haux, rots_cache, E1
        )
        println("Test grad psiks before rotate: $(2*dot(g, psiks))")
        println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
        rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
        println("Test grad psiks after rotate: $(2*dot(g, psiks))")
        println("Test grad Haux after rotate: $(dot(Haux, g_Haux))")
 
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
            @warn "Line minimization is not successful"
        end

        if Nspin == 2
            magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe integ magn = $magn")
        end
        ΔE = abs(E_new - E1)
        println("iterSD=$(iterSD) E_new = $(E_new), ΔE = $(ΔE)")
        println("Focc = ")
        display(Ham.electrons.Focc); println()
        println("ebands (w.r.t) Fermi energy = ")
        display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()

        found_zero_Focc = false
        for ikspin in 1:Nkspin
            s = sum(Ham.electrons.Focc[:,ikspin])
            if s < 1e-10
                @warn "Detected Focc for ikspin=$(ikspin) is zeros, break from iteration"
                found_zero_Focc = true
                #break # from ikspin loop
            end
        end
        if found_zero_Focc
            @warn "Some problems with Focc detected"
            #break
        end

        if ΔE < 1e-6
            println("Converged")
            break
        end 
        
        # New iterations, variables are updated in linmin_quad_v01
        E1 = E_new
    end

    @infiltrate

    return
end

