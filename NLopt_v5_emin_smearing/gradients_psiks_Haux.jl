
# for psiks
function calc_grad!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    g::BlochWavefunc,
    Hsub::Vector{Matrix{ComplexF64}}
)
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Focc = Ham.electrons.Focc

    #XXX Preallocate Hpsi ? Using Ngwx?
    # Ngwx = maximum(Ham.pw.gvecw.Ngw)

    for ispin in 1:Nspin, ik in 1:Nkpt
        # Don't forget to set current index for Hamiltonian
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        #
        Hpsi = op_H( Ham, psiks[ikspin] )
        Hsub[ikspin][:,:] = psiks[ikspin]' * Hpsi
        Hpsi[:,:] -= psiks[ikspin] * Hsub[ikspin]
        for ist in 1:Nstates
            g[ikspin][:,ist] .= Focc[ist,ikspin] .* Hpsi[:,ist]
        end
    end
    #
    return
end


# Gradient for Haux
# The real input is actually stored in Ham.electrons.ebands which is calculated
# from diagonalizing Haux
function calc_grad_Haux!(
    Ham, Hsub, g_Haux, Kg_Haux;
    κ=1.0
)

    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    fprime = zeros(Float64, Nstates)
    fprimeNum = zeros(Float64, Nstates)
    dmuNum = zeros(Float64, Nspin)
    dmuDen = zeros(Float64, Nspin)

    # These variables are not updated or calculated here
    # They are assumed to be calculated elsewhere
    ebands = Ham.electrons.ebands
    kT = Ham.electrons.kT
    E_fermi = Ham.electrons.E_fermi

    if Nspin == 1
        w_spin = 2.0
    else
        w_spin = 1.0 
    end
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        # accumulate with Nkpt? Add wk ?
        for ist in 1:Nstates
            fprime[ist] += smear_fermi_prime( ebands[ist,ikspin], E_fermi, kT )
            fprimeNum[ist] += fprime[ist] * ( real(Hsub[ikspin][ist,ist]) - ebands[ist,ikspin] )
        end
        # smear_fermi_prime might return NaN if E_fermi is not set properly
        dmuNum[ispin] += w_spin * sum(fprimeNum)
        dmuDen[ispin] += w_spin * sum(fprime)
    end

    dmuContrib = sum(dmuNum)/sum(dmuDen)
    dBzContrib = 0.0 # not used

    gradF0 = zeros(Float64, Nstates, Nstates)
    gradF = zeros(Float64, Nstates, Nstates)
    g_tmp = zeros(Float64, Nstates, Nstates)

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        gradF0[:,:] = real.(Hsub[ikspin]) - diagm( 0 => ebands[:,ikspin] )
        gradF[:,:] = copy(gradF0)
        for ist in 1:Nstates
            gradF[ist,ist] = gradF0[ist,ist] - dmuContrib # FIXME: not tested for spinpol
        end
        g_tmp[:,:] = grad_smear( smear_fermi, smear_fermi_prime, ebands[:,ikspin], E_fermi, kT, gradF )
        g_Haux[ikspin][:,:] = w_spin * 0.5 * (g_tmp' + g_tmp)
        Kg_Haux[ikspin][:,:] = -κ * gradF0[:,:] # preconditioning here?
    end

    return

end

#=
# Eq (24)
function calc_dFdmu(
    Hsub::Matrix{Float64},
    ebands,
    Focc_in,
    kT
)
    Nstates = size(ebands, 1)
    dFdmu = zeros(Float64, Nstates)
    if Nspin == 1
        Focc = 0.5*Focc_in
    else
        Focc = Focc_in
    end
    for ist in 1:Nstates
        dFdmu[ist] = (Hsub[ist,ist] - ebands[ist,1])*Focc[ist,1]*(1 - Focc[ist,1])
    end
    return dFdmu/kT
end


function offdiag_elements( Hsub, ebands, E_f::Float64, kT::Float64 )
    Nstates = size(evals, 1)
    mat = zeros(Nstates,Nstates)
    for j in 1:Nstates, i in 1:Nstates
        de = ebands[i] - ebands[j]
        if abs(de) > 1e-6
            mat[i,j] = Hsub[i,j] * ( smear_fermi(ebands[i], E_f, kT) - smear_fermi(ebands[j], E_f, kT) ) / de
        end
    end
    return mat
end
=#



function calc_grad_Lfunc_Haux!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Haux::Vector{Matrix{Float64}},
    g::BlochWavefunc,
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
        ebands[:,ikspin], Urot[:,:] = eigen(Hermitian(Haux[ikspin]))
        psiksU[ikspin] = psiks[ikspin]*Urot
    end

    update_from_ebands!( Ham, ebands )
    update_from_wavefunc!( Ham, psiksU )

    for ikspin in 1:Nkspin
        fill!(g[ikspin], 0.0 + 0.0*im)
        fill!(Hsub[ikspin], 0.0 + 0.0*im)
        fill!(g_Haux[ikspin], 0.0)
        fill!(Kg_Haux[ikspin], 0.0)
    end

    # Evaluate the gradients for psi
    calc_grad!(Ham, psiksU, g, Hsub) # don't forget to include Urot in psi
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
    calc_grad!(Ham, psi, g, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    return
end