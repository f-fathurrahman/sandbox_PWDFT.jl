#
# Various functions to update Hamiltonian
#

# Ham.electrons.ebands are already updated elsewhere
function update_from_ebands!(Ham)
    update_from_ebands!(Ham, Ham.electrons.ebands)
    return
end


# Input: ebands
# Modifies: Focc, E_fermi, mTS
# Ham.electrons.ebands are not modified here
# Also set kT (hardcoded)
function update_from_ebands!(Ham, ebands)

    # NOTE: ebands are assumed to be updated outside this function

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    Nelectrons = Ham.electrons.Nelectrons
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    kT = Ham.electrons.kT
    @assert kT > 1e-10

    E_fermi, mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT,
        Nkpt, wk
    )
    # Set some output
    Ham.electrons.E_fermi = E_fermi
    Ham.energies.mTS = mTS

    return
end


# Input: psiks
# Modifies: Ham.rhoe, potentials
function update_from_wavefunc!(Ham, psiks)    
    # Compute electron density from psiks
    # Use Ham.rhoe
    calc_rhoe!(Ham, psiks, Ham.rhoe)
    # Update the potentials
    update_from_rhoe!(Ham, psiks, Ham.rhoe)
    # XXX: update_from_rhoe! will not overwrite update Ham.rhoe
    return
end

function get_diag_Haux_from_ebands( Ham )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
    end
    return Haux
end

# Modifies Ham.electrons.ebands
function transform_psiks_Haux_update_ebands!(Ham, psiks, Haux)
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    ebands = Ham.electrons.ebands
    Urot = zeros(ComplexF64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        ebands[:,ikspin], Urot[:,:] = eigen(Hermitian(Haux[ikspin]))
        psiks[ikspin] *= Urot # rotate psiks
        Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
    end
    return
end

# Haux (in diagonal form) is stored in Ham.electrons.ebands
# psiks is already orthonormalized and rotated according to Urot
# that makes Haux diagonal
# 
# The following must be 
# update_from_ebands!( Ham, ebands )
# update_from_wavefunc!( Ham, psiks )
#
function calc_Lfunc(
    Ham::Hamiltonian,
    psiks::BlochWavefunc
)
    calc_energies!(Ham, psiks)
    # get entropy
    # Ham.electrons.mTS is computed in update_from_ebands!
    Ham.energies.mTS = Ham.electrons.mTS
    #
    return sum(Ham.energies)
end

# The inputs are:
# - wavefunction psi, and
# - auxiliary Hamiltonian in diagonal form, stored as matrix with size (Nstates,Nspin)
#
# Some fields of Ham will be modified
function calc_Lfunc_ebands!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    ebands::Matrix{Float64} # (Nstates,Nkspin)
)
    update_from_ebands!( Ham, ebands )
    update_from_wavefunc!( Ham, psiks )
    #
    calc_energies!(Ham, psiks)
    # get entropy
    # Ham.electrons.mTS is computed in update_from_ebands!
    Ham.energies.mTS = Ham.electrons.mTS
    
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
    Haux::Vector{Matrix{Float64}}
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

    return calc_Lfunc_ebands!(Ham, psiksU, ebands)
end

