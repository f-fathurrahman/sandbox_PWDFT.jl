using Revise, Infiltrator

using LinearAlgebra
using Printf
using Random

using PWDFT

includet("smearing.jl")
includet("occupations.jl")
includet("Lfunc.jl")
includet("gradients_psiks_Haux.jl")
includet("prepare_Ham_various.jl")

# Ham.electrons.ebands and Ham.rhoe are modified
# psiks and Haux are not modified
function try_step!(α::Float64, Ham, psiks, Haux, d, d_Haux)
    Nkspin = length(psiks)
    # Step
    psiks_new = psiks + α*d
    Haux_new = Haux + α*d_Haux
    for ikspin in 1:Nkspin
        ortho_sqrt!(Ham, psiks_new[ikspin])
    end
    transform_psiks_Haux_update_ebands!( Ham, psiks_new, Haux_new )
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks_new )
    #
    E_try = calc_Lfunc( Ham, psiks_new )
    return E_try
end


function main_sd_01(Ham; NiterMax=100, α=0.1)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Rhoe = Ham.rhoe
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

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

    # Prepare Haux (random numbers)
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    # eigenvalues of Hsub
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm(0 => eigvals(Hsub[ikspin]))
    end

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

    if Nspin == 2
        magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
        integRhoe = sum(Rhoe)*dVol
        println("INITIAL INPUT: integRhoe = $integRhoe integ magn = $magn")
    end

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux )
    #
    # Update Hamiltonian before evaluating free energy
    update_from_ebands!( Ham )
    #update_from_wavefunc!( Ham, psiks )
    E1 = calc_Lfunc( Ham, psiks )
    @info "E1 from Lfunc_Haux = $(E1)"

    #calc_grad_Lfunc_Haux!( Ham, psiks, Haux, g, Hsub, g_Haux, Kg_Haux )

    calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

    println("Test grad psiks: $(2*dot(g, psiks))")
    println("Test grad Haux: $(dot(Haux, g_Haux))")

    if Nspin == 2
        magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
        integRhoe = sum(Rhoe)*dVol
        println("INITIAL after evaluate: integRhoe = $integRhoe integ magn = $magn")
    end

    for iterSD in 1:NiterMax

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin] = -Kg[ikspin]
            d_Haux[ikspin] = -Kg_Haux[ikspin]
        end
        constrain_search_dir!(d, psiks)

        gd = 2*real(dot(g,d)) + real(dot(g_Haux, d_Haux))
        @info "gd = $(gd)"
        if gd > 0
            error("Bad step direction")
        end

        α_t = 5.0
        is_linmin_success = false
        for itry in 1:10
            E_t = try_step!(α_t, Ham, psiks, Haux, d, d_Haux)
            println("itry=$(itry), using α_t=$(α_t) E_t=$(E_t)")
            if E_t < E1
                println("α_t is accepted")
                is_linmin_success = true
                break
            end
            α_t = α_t * 0.5 # reduce
        end

        if !is_linmin_success
            @warn "LineMin is not successful, break from iteration"
            break
        end

        # Step
        psiks .+= α_t*d
        Haux .+= α_t*d_Haux

        for ikspin in 1:Nkspin
            ortho_sqrt!(Ham, psiks[ikspin])
        end
        transform_psiks_Haux_update_ebands!( Ham, psiks, Haux )
        update_from_ebands!( Ham )
        update_from_wavefunc!( Ham, psiks )
        #
        E2 = calc_Lfunc( Ham, psiks )
        #
        calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
        println("Test grad psiks: $(2*dot(g, psiks))")
        println("Test grad Haux: $(dot(Haux, g_Haux))")

        #
        ΔE = E2 - E1
        println("iterSD=$(iterSD) E=$(E2) abs(ΔE)=$(abs(ΔE)) E_fermi=$(Ham.electrons.E_fermi)")
        #
        println("Focc = ")
        display(Ham.electrons.Focc); println()
        println("ebands (w.r.t) Fermi energy = ")
        display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()
        #
        if Nspin == 2
            magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe integ magn = $magn")
        end

        #
        if ΔE > 0
            @warn "Energy is increasing"
            break
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
