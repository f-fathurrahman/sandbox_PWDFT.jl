using LinearAlgebra
using Printf
using Infiltrator
using Serialization

using PWDFT

include("smearing.jl")
include("occupations.jl")
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

    # This is only for computing diagonal Haux
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

    # Compute this once
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    E1 = calc_Lfunc_Haux!( Ham, psiks, Haux )
    @info "E1 from Lfunc_Haux = $(E1)"

    calc_grad_Lfunc_Haux!( Ham, psiks, Haux, g, Hsub, g_Haux, Kg_Haux )

    # Compute Haux again
    Hsub_after = deserialize("Hsub_after.jldat")
    calc_grad_Haux!( Ham, Hsub_after, g_Haux, Kg_Haux )


    @infiltrate

end

main()


