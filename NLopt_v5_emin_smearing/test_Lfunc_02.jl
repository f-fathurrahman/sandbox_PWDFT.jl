using LinearAlgebra
using Printf

using PWDFT

include("smearing.jl")
include("occupations.jl")
include("update_Hamiltonian.jl")
include("Lfunc.jl")
include("utilities_emin_smearing.jl")

function main()

    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");

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


    Haux = Vector{Matrix{Float64}}(undef, Nkspin)
    Urot = zeros(Float64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(Float64, Nstates, Nstates)
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    end

    F1 = calc_Lfunc_Haux!( Ham, psiks, Haux )
    @info "F1 from Lfunc_Haux = $(F1)"

    Urot = zeros(ComplexF64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        Urot[:,:] = random_unitary_matrix(Nstates)
        psiks[ikspin][:,:] = psiks[ikspin] * Urot
        Haux[ikspin][:,:] = Urot' * Haux[ikspin] * Urot
    end
    F2 = calc_Lfunc_Haux!( Ham, psiks, Haux )
    @info "F2 from Lfunc_ebands = $(F2)"
end

main()
