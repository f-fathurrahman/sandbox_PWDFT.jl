using LinearAlgebra
using Printf
using Infiltrator
using Random

using PWDFT

include("smearing.jl")
include("occupations.jl")
include("update_Hamiltonian.jl")
include("Lfunc.jl")
include("gradients_psiks_Haux.jl")
include("utilities_emin_smearing.jl")

function main()

    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");

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

    # Prepare Haux
    Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(Float64, Nstates, Nstates)
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    end

    g = zeros_BlochWavefunc(Ham)
    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    g_Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        g_Haux[ikspin] = zeros(Float64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(Float64, Nstates, Nstates)
    end

    # Compute this once
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    # Make Haux diagonal
    transform_psiks_Haux!( Ham, psiks, Haux )

    E1 = calc_Lfunc_Haux!( Ham, psiks, Haux )
    @info "E1 from Lfunc_Haux = $(E1)"

    calc_grad_Lfunc_Haux!( Ham, psiks, Haux, g, Hsub, g_Haux, Kg_Haux )

    Δ = 1e-5
    Δ_Haux = 1e-5
    dW = Δ*g
    dW_Haux = Δ_Haux*g_Haux
    psiks_new = psiks + dW
    Haux_new = Haux + dW_Haux
    for ikspin in 1:Nkspin
        ortho_sqrt!(Ham, psiks_new[ikspin])
    end
    E2 = calc_Lfunc_Haux!( Ham, psiks_new, Haux_new )
    @info "E2 from Lfunc_Haux = $(E2)"

    ΔE = abs(E1 - E2)
    dE_psiks = 2*real(dot(g, dW))
    dE_Haux = dot(g_Haux, dW_Haux)
    println("dE_psiks = ", dE_psiks)
    println("dE_Haux = ", dE_Haux)
    println("ratio abs(ΔE)/dE = ", ΔE/(dE_psiks + dE_Haux) )

    @infiltrate

end



main()
