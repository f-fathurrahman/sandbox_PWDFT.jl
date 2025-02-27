using LinearAlgebra
using Printf
using Serialization
using Infiltrator

using PWDFT

include("smearing.jl")
include("occupations.jl")
include("gradients_psiks_Haux.jl")
include("Lfunc.jl")
include("utilities_emin_smearing.jl")

function main()

    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");

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

    psiks = deserialize("DATA_jdftx_Fe_v2/psiks.jldat")
    Focc_r = deserialize("DATA_jdftx_Fe_v2/Focc.jldat")

    Hsub = deserialize("DATA_jdftx_Fe_v2/Hsub.jldat")
    Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = zeros(Float64, Nstates, Nstates)
        Haux[ikspin][:,:] = diagm( 0 => eigvals(Hermitian(Hsub[ikspin])) )
    end

    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end

    # Compute this once
    Ham.energies.NN = calc_E_NN(Ham.atoms)

    F1 = calc_Lfunc_Haux!( Ham, psiks, Haux )
    @info "F1 from Lfunc_Haux = $(F1)"

    # Also compute 
    calc_grad_psiks!(Ham, psiks, g, Hsub)
    my_Kprec!(Ham, g, Kg) # not really needed
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

    g_Haux_r = deserialize("DATA_jdftx_Fe_v2/g_Haux.jldat")
    Kg_Haux_r = deserialize("DATA_jdftx_Fe_v2/Kg_Haux.jldat")

    @infiltrate

end

main()
