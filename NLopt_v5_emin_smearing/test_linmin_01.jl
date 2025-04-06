using LinearAlgebra
using Printf
using Infiltrator
using Random

using PWDFT

include("smearing.jl")
include("occupations.jl")
include("Lfunc.jl")
include("gradients_psiks_Haux.jl")
include("utilities_emin_smearing.jl")

function main()

    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    Random.seed!(1234)

    # This will take into account whether the overlap operator is needed or not
    psiks = rand_BlochWavefunc(Ham)

    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
        Ham.electrons.kT = kT
    end
    # No need to set kT here, it is already default to 0

    if pwinput.nspin == 2
        starting_magnetization = pwinput.starting_magnetization
    else
        starting_magnetization = nothing
    end

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk
    Nelectrons = Ham.electrons.Nelectrons

    # Prepare Haux (random numbers)
    Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(Float64, Nstates, Nstates)
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    end

    # Gradients, subspace Hamiltonian
    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    d = zeros_BlochWavefunc(Ham)
    #
    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    g_Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    d_Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        g_Haux[ikspin] = zeros(Float64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(Float64, Nstates, Nstates)
        d_Haux[ikspin] = zeros(Float64, Nstates, Nstates)
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
    println("Initial energy = $(E1)")
    # Gradients
    calc_grad_psiks!(Ham, psiks, g, Hsub)
    my_Kprec!(Ham, g, Kg)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    #
    println("Test grad psiks: $(2*dot(g, psiks))")
    println("Test grad Haux: $(dot(Haux, g_Haux))")


    for iterSD in 1:30

        println("\nBegin iterSD = $(iterSD)")

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin][:,:] = -Kg[ikspin][:,:]
            d_Haux[ikspin][:,:] = -Kg_Haux[ikspin][:,:]
        end
        constrain_search_dir!(Ham, d, psiks)

        α = 1.0 #linmin_quad!(Ham, psiks, Haux, g, g_Haux, d, d_Haux, E1)
        println("Using α=$(α) from linmin")

        # Step
        psiks .+= α*d
        Haux .+= α*d_Haux

        for ikspin in 1:Nkspin
            ortho_sqrt!(Ham, psiks[ikspin])
        end
        transform_psiks_Haux_update_ebands!( Ham, psiks, Haux )
        update_from_ebands!( Ham )
        update_from_wavefunc!( Ham, psiks )
        #
        E2 = calc_Lfunc( Ham, psiks )
        #
        calc_grad_psiks!(Ham, psiks, g, Hsub)
        my_Kprec!(Ham, g, Kg)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

        println("Test grad psiks: $(2*dot(g, psiks))")
        println("Test grad Haux: $(dot(Haux, g_Haux))")

        #
        ΔE = E2 - E1
        println("iterSD=$(iterSD) E=$(E2) abs(ΔE)=$(abs(ΔE)) E_fermi=$(Ham.electrons.E_fermi)")
        #
        if ΔE > 0
            @warn "Energy is increasing"
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


function linmin_quad!(Ham, psiks, Haux, g, g_Haux, d, d_Haux, E_old)

    gd = 2*real(dot(g,d)) + dot(g_Haux, d_Haux)
    println("gd = $(gd)")
    if gd > 0
        @infiltrate
        error("Bad step direction")
    end

    NtryMax = 5

    α_t = 1.0
    α_safe =  1e-5 # safe step size
    α = α_safe
    IS_SUCCESS = false
    for itry in 1:NtryMax
        println("--- Begin itry linmin = $(itry) using α_t=$(α_t)")
        E_t = try_step!(α_t, Ham, psiks, Haux, d, d_Haux)
        c = ( E_t - (E_old + α_t*gd) ) / α_t^2
        α = -gd/(2*c)
        println("Find α = $(α)")
        if α < 0
            @info "Wrong curvature, α is negative"
            α_t *= 2.0
            println("Should continue to next iter")
            continue
        end
        #
        E_t2 = try_step!(α, Ham, psiks, Haux, d, d_Haux)
        println("Trial energy 2: E_t2 = $(E_t2)")
        if E_t2 > E_old
            @warn "Energy is not decreasing"
            α_t *= 0.1
        else
            println("Accept α=$(α)")
            IS_SUCCESS = true
            break
        end
    end

    if IS_SUCCESS
        return α
    else
        return α_safe
    end

end


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

main()
