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
    Haux::Vector{Matrix{ComplexF64}}
)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt * Nspin
    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands

    psiksU = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Urot = zeros(ComplexF64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        ebands[:,ikspin], Urot[:,:] = eigen(Hermitian(Haux[ikspin]))
        psiksU[ikspin] = psiks[ikspin]*Urot
    end

    return calc_Lfunc_ebands!(Ham, psiksU, ebands)
end


function calc_grad_Lfunc_Haux!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Haux::Vector{Matrix{Float64}},
    g::BlochWavefunc, Kg::BlochWavefunc,
    Hsub,
    g_Haux,
    Kg_Haux
)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt * Nspin
    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands

    # Diagonalize Haux and store the results in ebands
    # Also rotate psiks accordingly
    psiksU = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Urot = zeros(Float64, Nstates, Nstates)
    for ikspin in 1:Nkspin
        # diagonal Haux is stored as ebands
        ebands[:,ikspin], Urot[:,:] = eigen(Hermitian(Haux[ikspin]))
        # Also rotate psiks
        psiksU[ikspin] = psiks[ikspin]*Urot
    end

    update_from_ebands!( Ham, ebands )
    update_from_wavefunc!( Ham, psiksU )

    # Evaluate the gradients for psi
    calc_grad_psiks!(Ham, psiksU, g, Kg, Hsub) # don't forget to include Urot in psi
    # pass Hsub
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

    return
end

function calc_grad_Lfunc_ebands!(
    Ham::Hamiltonian,
    psi, # (Nbasis,Nstates)
    ebands, # (Nstates,1)
    g,
    Hsub,
    g_Haux,
    Kg_Haux
)

    @assert size(ebands,2) == 1

    Haux = diagm(0 => ebands[:,1])

    update_from_ebands!(Ham, ebands)
    update_from_wavefunc!(Ham, psi)

    fill!(g, 0.0)
    fill!(Hsub, 0.0)
    fill!(g_Haux, 0.0)
    fill!(Kg_Haux, 0.0)

    # Evaluate the gradient for psi
    calc_grad_psiks!(Ham, psi, g, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    return
end

