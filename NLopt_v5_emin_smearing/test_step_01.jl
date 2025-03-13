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


mutable struct RotationsCache
    rotPrev::Vector{Matrix{ComplexF64}}
    rotPrevC::Vector{Matrix{ComplexF64}}
    rotPrevCinv::Vector{Matrix{ComplexF64}}
    Urot::Vector{Matrix{ComplexF64}}
    UrotC::Vector{Matrix{ComplexF64}}
end

function RotationsCache(Nkspin, Nstates)
    rotPrev = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Urot = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    UrotC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        Urot[ikspin] = Matrix(1.0*I(Nstates))
        UrotC[ikspin] = Matrix(1.0*I(Nstates))
    end
    return RotationsCache(
        rotPrev,
        rotPrevC,
        rotPrevCinv,
        Urot,
        UrotC,
    )
end


function rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
    Nkspin = length(g)
    rotPrevCinv = rots_cache.rotPrevCinv
    rotPrev = rots_cache.rotPrev
    for ikspin in 1:Nkspin
        g[ikspin][:,:] = g[ikspin][:,:] * rotPrevCinv[ikspin]
        Kg[ikspin][:,:] = Kg[ikspin][:,:] * rotPrevCinv[ikspin]
        g_Haux[ikspin][:,:] = rotPrev[ikspin] * g_Haux[ikspin][:,:] * rotPrev[ikspin]'
        Kg_Haux[ikspin][:,:] = rotPrev[ikspin] * Kg_Haux[ikspin][:,:] * rotPrev[ikspin]'
    end
    return
end

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
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
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

    rots_cache = RotationsCache(Nkspin, Nstates)

    Haux_orig = copy(Haux)

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache,
        do_ortho_psi=false, overwrite_Haux=true
    )
    # rots_cache is needed because Haux is not yet diagonal
    # do_ortho_psi=false because psiks is already orthonormal, UrotC is identity
    #
    # Update Hamiltonian before evaluating free energy
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks )
    E1 = calc_Lfunc( Ham, psiks )
    @info "E1 from Lfunc_Haux = $(E1)"

    # Calculate gradients
    calc_grad_psiks!(Ham, psiks, g, Hsub)
    my_Kprec!(Ham, g, Kg)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    #
    println("Test grad psiks before rotate: $(2*dot(g, psiks))")
    println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
    #
    #
    # Gradients are also rotated
    # (XXX: do we absolutely need this?)
    # XXX: If this is used, then the gradients are w.r.t variables (psiks and Haux) before transformation
    # XXX: I think we don't need this
    #=
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
    #
    println("Test grad psiks after rotate: $(2*dot(g, psiks))")
    println("Test grad Haux after rotate: $(dot(Haux, g_Haux))")
    #
    println("Test grad Haux_orig after rotate: $(dot(Haux_orig, g_Haux))")
    =#
    # need to save Haux_eigs ?

    #@infiltrate

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

    α = 1.0
    # Step forward (positive α)
    E2 = do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)
    println("E2 = ", E2)
    # E2 should be lower (if not overshooting linemin)
    if E1 < E2
        @warn "E2 is larger"
    end

    # Step backward (negative α)
    E3 = do_step_psiks_Haux!(-α, Ham, psiks, Haux, d, d_Haux, rots_cache)
    println("E3 = ", E3)
    println("E3 should be the same as E1 = ", E1)

    # TODO: need to compare variables psiks and Haux after step backward and original

    #@infiltrate

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


# Ham.electrons.ebands and Ham.rhoe are modified
# psiks and Haux are modified
function do_step_psiks_Haux!(α::Float64, Ham, psiks, Haux, d, d_Haux, rots_cache)

    Nkspin = length(psiks)

    ebands = Ham.electrons.ebands

    UrotC = rots_cache.UrotC
    Urot = rots_cache.Urot
    rotPrev = rots_cache.rotPrev
    rotPrevC = rots_cache.rotPrevC
    rotPrevCinv = rots_cache.rotPrevCinv

    # Step
    for ikspin in 1:Nkspin
        psiks[ikspin] += α * d[ikspin] * rotPrevC[ikspin]
        Haux[ikspin]  += α * rotPrev[ikspin]' * d_Haux[ikspin] * rotPrev[ikspin]
    end

    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache, do_ortho_psi=true )

    # Update Hamiltonian terms
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks )
    # Now, we are ready to evaluate
    E_try = calc_Lfunc( Ham, psiks )
    return E_try
end


main()
