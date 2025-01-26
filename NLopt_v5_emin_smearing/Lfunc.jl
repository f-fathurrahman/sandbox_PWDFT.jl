include("update_Hamiltonian.jl")

# The inputs are:
# - wavefunction psi, and
# - auxiliary Hamiltonian in diagonal form, stored as matrix with size (Nstates,Nspin)
#
# Some fields of Ham will be modified
function calc_Lfunc_ebands!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    ebands::Matrix{Float64}, # (Nstates,Nkspin)
    kT::Float64
)
    E_fermi, mTS = update_from_ebands!( Ham, ebands, kT )
    #@info "E_fermi = $(E_fermi)"
    update_from_wavefunc!( Ham, psiks )
    #
    calc_energies!(Ham, psiks)
    Ham.energies.mTS = mTS
    
    return sum(Ham.energies)
end


# The inputs are:
# - wavefunction psi, and
# - auxiliary Hamiltonian Haux. No support for spin polarization for now.
#
# Some fields of Ham will be modified
#
# psi and Haux must be transformed simultaneously by using some unitary matrix.
# The transformation chosen such that Haux transformed to diagonal form using
# eigendecomposition.
#
# psi and Haux will not be modified in place upon calling this function.
function calc_Lfunc_Haux!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Haux::Vector{Matrix{Float64}},
    kT
)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt * Nspin
    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands

    psiksU = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Urot = zeros(Float64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        ebands[:,ikspin], Urot[:,:] = eigen(Hermitian(Haux[ikspin]))
        psiksU[ikspin] = psiks[ikspin]*Urot
    end

    return calc_Lfunc_ebands!(Ham, psiksU, ebands, kT)
end

