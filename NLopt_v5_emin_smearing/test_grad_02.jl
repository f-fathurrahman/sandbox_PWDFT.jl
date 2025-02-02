using LinearAlgebra
using Printf
using Infiltrator

using PWDFT

include("smearing.jl")
include("occupations.jl")
include("update_Hamiltonian.jl")
include("Lfunc.jl")
include("gradients_psiks_Haux.jl")

function main()

    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");
    # psiks and Haux will be read from files

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
    CellVolume = Ham.pw.CellVolume

    psiks = deserialize("psiks.jldat")
    # Renormalize
    for ikspin in 1:Nkspin
        psiks[ikspin] .*= sqrt(CellVolume)
    end

    # Prepare Haux
    Hsub = deserialize("Hsub.jldat")
    Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = zeros(Float64, Nstates, Nstates)
        Haux[ikspin][:,:] = diagm( 0 => eigvals(Hermitian(Hsub[ikspin])) )
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

    E1 = calc_Lfunc_Haux!( Ham, psiks, Haux )
    @info "E1 from Lfunc_Haux = $(E1)"

    calc_grad_Lfunc_Haux!( Ham, psiks, Haux, g, Hsub, g_Haux, Kg_Haux )

    #=
    Δ = 1e-5
    psiks_new = psiks + Δ*g;
    Haux_new = Haux #+ Δ*g_Haux;
    for ikspin in 1:Nkspin
        ortho_sqrt!(Ham, psiks_new[ikspin])
    end
    E2 = calc_Lfunc_Haux!( Ham, psiks_new, Haux_new )
    @info "E2 from Lfunc_Haux = $(E2)"
    =#

    @infiltrate
    @exfiltrate

end

main()
