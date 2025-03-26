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

function main_sd_02(Ham; NiterMax=100)

    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands

    # Prepare Haux (random numbers), diagonal
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm(0 => sort(randn(Float64, Nstates)))
    end

    # Gradients, subspace Hamiltonian
    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    d = zeros_BlochWavefunc(Ham)
    #
    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    d_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        d_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux )
    #
    # Update Hamiltonian before evaluating free energy
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks )
    E1 = calc_Lfunc( Ham, psiks )
    println("E1 from Lfunc_Haux = $(E1)")

    #calc_grad_Lfunc_Haux!( Ham, psiks, Haux, g, Hsub, g_Haux, Kg_Haux )

    calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

    println("Test grad psiks: $(2*dot(g, psiks))")
    println("Test grad Haux: $(dot(Haux, g_Haux))")

    rotPrev = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    #
    # Preallocated
    Urot = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    UrotC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        Urot[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        UrotC[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end

    α = 0.1
    Rhoe = Ham.rhoe
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    for iterSD in 1:NiterMax

        println("\nBegin SD iteration=$(iterSD)")

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin][:,:] = -Kg[ikspin][:,:]
            d_Haux[ikspin][:,:] = -Kg_Haux[ikspin][:,:]
        end
        constrain_search_dir!(d, psiks)

        gd = 2*real(dot(g,d)) + real(dot(g_Haux, d_Haux))
        @info "gd = $(gd)"
        if gd > 0
            error("Bad step direction")
        end

        # Step
        for ikspin in 1:Nkspin
            psiks[ikspin] += α * d[ikspin] * rotPrevC[ikspin]
            Haux[ikspin]  += α * rotPrev[ikspin]' * d_Haux[ikspin] * rotPrev[ikspin]
        end

        # Transform
        for ikspin in 1:Nkspin
            ebands[:,ikspin], Urot[ikspin][:,:] = eigen(Hermitian(Haux[ikspin]))
            Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
            # wavefunc
            UrotC[ikspin][:,:] = inv(sqrt(psiks[ikspin]' * psiks[ikspin]))
            UrotC[ikspin][:,:] *= Urot[ikspin] # extra rotation
            psiks[ikspin][:,:] = psiks[ikspin]*UrotC[ikspin]
        end

        for ikspin in 1:Nkspin
            rotPrev[ikspin] = rotPrev[ikspin] * Urot[ikspin]
            rotPrevC[ikspin] = rotPrevC[ikspin] * UrotC[ikspin]
            rotPrevCinv[ikspin] = inv(UrotC[ikspin]) * rotPrevCinv[ikspin]
        end

        update_from_ebands!( Ham )
        update_from_wavefunc!( Ham, psiks )
        #
        E2 = calc_Lfunc( Ham, psiks )
        #
        calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

        println("Test grad psiks before rotate: $(2*dot(g, psiks))")
        println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
        for ikspin in 1:Nkspin
            g[ikspin][:,:] = g[ikspin][:,:] * rotPrevCinv[ikspin]
            Kg[ikspin][:,:] = Kg[ikspin][:,:] * rotPrevCinv[ikspin]
            g_Haux[ikspin][:,:] = rotPrev[ikspin] * g_Haux[ikspin][:,:] * rotPrev[ikspin]'
            Kg_Haux[ikspin][:,:] = rotPrev[ikspin] * Kg_Haux[ikspin][:,:] * rotPrev[ikspin]'
        end
        println("Test grad psiks AFTER rotate: $(2*dot(g, psiks))")
        println("Test grad Haux AFTER rotate: $(dot(Haux, g_Haux))")

        #
        ΔE = E2 - E1
        @info "iterSD=$(iterSD) E=$(E2) abs(ΔE)=$(abs(ΔE)) E_fermi=$(Ham.electrons.E_fermi)"
        if Nspin == 2
            magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe integ magn = $magn")
        end
        println("Focc = ")
        display(Ham.electrons.Focc); println
        println("ebands (w.r.t Fermi) = ")
        display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println

        #
        if ΔE > 0
            error("Energy is increasing")
        end
        #
        if abs(ΔE) < 1e-6
            @info "Converged"
            break
        end

        E1 = E2

    end

    @infiltrate

end
