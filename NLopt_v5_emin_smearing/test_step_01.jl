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
    @info "E1 from Lfunc_Haux = $(E1)"

    #calc_grad_Lfunc_Haux!( Ham, psiks, Haux, g, Hsub, g_Haux, Kg_Haux )

    calc_grad_psiks!(Ham, psiks, g, Hsub)
    my_Kprec!(Ham, g, Kg)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

    println("Test grad psiks: $(2*dot(g, psiks))")
    println("Test grad Haux: $(dot(Haux, g_Haux))")


    # Set direction
    for ikspin in 1:Nkspin
        d[ikspin][:,:] = -Kg[ikspin][:,:]
        d_Haux[ikspin][:,:] = -Kg_Haux[ikspin][:,:]
    end
    constrain_search_dir!(d, psiks)

    gd = 2*real(dot(g,d)) + dot(g_Haux, d_Haux)
    @info "gd = $(gd)"
    if gd > 0
        error("Bad step direction")
    end

    α = 1.0
    E_old = E1
    E_new = try_step!(α, Ham, psiks, Haux, d, d_Haux)
    println("E_new = ", E_new)
    if E_old < E_new
        @warn "E_new is larger"
    end

    E_old = E_new
    E_new = try_step!(α, Ham, psiks, Haux, d, d_Haux)
    println("E_new = ", E_new)
    if E_old < E_new
        @warn "E_new is larger"
    end

    @infiltrate

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
