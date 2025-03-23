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


function main()

    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");
    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms);
    
    # We need to set some parameters manually:
    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
        Ham.electrons.kT = kT
    end

    if pwinput.nspin == 2
        starting_magnetization = pwinput.starting_magnetization
    else
        starting_magnetization = nothing
    end

    # %%
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates;

    # %% [markdown]
    # Initialize electronic variables: `psiks` and `Haux`:

    # %%
    Random.seed!(1234)
    # This will take into account whether the overlap operator is needed or not
    psiks = rand_BlochWavefunc(Ham);

    # Prepare Haux (random numbers)
    #
    # For Haux, choose between generic symmetric Haux:
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    end

    #=
    # or diagonal Haux:
    # Prepare Haux (random numbers)
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

    rots_cache = RotationsCache(Nkspin, Nstates);
    Haux_orig = copy(Haux);

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache,
        do_ortho_psi=false, overwrite_Haux=true
    )

    # Update Hamiltonian, compute energy and gradients at current psiks and Haux:

    # %%
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

    println("Test grad psiks before rotate: $(2*dot(g, psiks))")
    println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
    # %%
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
    # %%
    println("Test grad psiks after rotate: $(2*dot(g, psiks))")
    println("Test grad Haux after rotate: $(dot(Haux, g_Haux))")
    # %%
    println("Test grad Haux orig after rotate: $(dot(Haux_orig, g_Haux))")

    α_t_start = 1.0
    α_t_min = 1e-10
    α_t = α_t_start

    Rhoe = Ham.rhoe
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    for iterSD in 1:500

        println("\nStart iterSD = ", iterSD)

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin][:,:] = -Kg[ikspin][:,:]
            d_Haux[ikspin][:,:] = -Kg_Haux[ikspin][:,:]
        end
        constrain_search_dir!(d, psiks)

        #
        # Do line minimization:
        E_new, is_success, α = linmin_quad_v01!(
            α_t,
            Ham, psiks, Haux, Hsub, g, g_Haux, Kg, Kg_Haux, d, d_Haux, rots_cache, E1
        )
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
            @warn "Line minimization is not successful"
        end

        if Nspin == 2
            magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe integ magn = $magn")
        end
        ΔE = abs(E_new - E1)
        println("iterSD=$(iterSD) E_new = $(E_new), ΔE = $(ΔE)")

        if ΔE < 1e-6
            println("Converged")
            break
        end 
        
        # New iterations, variables are updated in linmin_quad_v01
        E1 = E_new
    end

end

