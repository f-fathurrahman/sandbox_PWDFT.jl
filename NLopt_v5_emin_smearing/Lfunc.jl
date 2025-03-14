mutable struct RotationsCache
    rotPrev::Vector{Matrix{ComplexF64}}
    rotPrevC::Vector{Matrix{ComplexF64}}
    rotPrevCinv::Vector{Matrix{ComplexF64}}
    Urot::Vector{Matrix{ComplexF64}}
    UrotC::Vector{Matrix{ComplexF64}}
end

function RotationsCache(Nkspin, Nstates)
    rotPrev = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Urot = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    UrotC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        Urot[ikspin] = Matrix(1.0*I(Nstates))
        UrotC[ikspin] = Matrix(1.0*I(Nstates))
    end
    return RotationsCache(
        rotPrev,
        rotPrevC,
        rotPrevCinv,
        Urot,
        UrotC,
    )
end


function rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
    Nkspin = length(g)
    rotPrevCinv = rots_cache.rotPrevCinv
    rotPrev = rots_cache.rotPrev
    for ikspin in 1:Nkspin
        g[ikspin][:,:] = g[ikspin][:,:] * rotPrevCinv[ikspin]
        Kg[ikspin][:,:] = Kg[ikspin][:,:] * rotPrevCinv[ikspin]
        g_Haux[ikspin][:,:] = rotPrev[ikspin] * g_Haux[ikspin][:,:] * rotPrev[ikspin]'
        Kg_Haux[ikspin][:,:] = rotPrev[ikspin] * Kg_Haux[ikspin][:,:] * rotPrev[ikspin]'
    end
    return
end



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



# Modifies Ham.electrons.ebands, save rotations in rots_cache
# psiks is assumed to be
function transform_psiks_Haux_update_ebands!(
    Ham, psiks, Haux, rots_cache;
    overwrite_Haux=true,
    do_ortho_psi=true
)
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    ebands = Ham.electrons.ebands
    #
    Urot = rots_cache.Urot
    UrotC = rots_cache.UrotC
    rotPrev = rots_cache.rotPrev
    rotPrevC = rots_cache.rotPrevC
    rotPrevCinv = rots_cache.rotPrevCinv
    #
    for ikspin in 1:Nkspin
        ebands[:,ikspin], Urot[ikspin] = eigen(Hermitian(Haux[ikspin]))
        if overwrite_Haux
            Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
        end
        #
        if do_ortho_psi
            UrotC[ikspin][:,:] = inv(sqrt(psiks[ikspin]' * psiks[ikspin]))
        end
        UrotC[ikspin][:,:] *= Urot[ikspin] # extra rotation
        psiks[ikspin][:,:] = psiks[ikspin]*UrotC[ikspin]
    end
    #
    for ikspin in 1:Nkspin
        # Save previous
        rotPrev[ikspin] *= Urot[ikspin]
        rotPrevC[ikspin] = rotPrevC[ikspin] * UrotC[ikspin]
        rotPrevCinv[ikspin] = inv(UrotC[ikspin]) * rotPrevCinv[ikspin]
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


function do_step_psiks_Haux!(α::Float64, Ham, psiks, Haux, d, d_Haux, rots_cache)

    Nkspin = length(psiks)

    ebands = Ham.electrons.ebands

    UrotC = rots_cache.UrotC
    Urot = rots_cache.Urot
    rotPrev = rots_cache.rotPrev
    rotPrevC = rots_cache.rotPrevC
    rotPrevCinv = rots_cache.rotPrevCinv

    # Step
    for ikspin in 1:Nkspin
        psiks[ikspin] += α * d[ikspin] * rotPrevC[ikspin]
        Haux[ikspin]  += α * rotPrev[ikspin]' * d_Haux[ikspin] * rotPrev[ikspin]
    end

    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache, do_ortho_psi=true )
    return
end



function do_compute_energy(Ham, psiks)
    # Update Hamiltonian terms
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks )
    # Now, we are ready to evaluate
    E = calc_Lfunc( Ham, psiks )
    return E
end



function linmin_quad_v01!(Ham, psiks, Haux, g, g_Haux, d, d_Haux, E_old)

    gd = 2*real(dot(g,d)) + real(dot(g_Haux, d_Haux))
    println("gd = $(gd)")
    if gd > 0
        @infiltrate
        error("Bad step direction")
    end

    NtryMax = 5

    α_t = 1.0
    α_safe =  1e-5 # safe step size
    α = α_safe
    IS_SUCCESS = false
    for itry in 1:NtryMax
        println("--- Begin itry linmin = $(itry) using α_t=$(α_t)")
        do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)
        E_t = do_compute_energy(Ham, psiks)
        c = ( E_t - (E_old + α_t*gd) ) / α_t^2
        α = -gd/(2*c)
        println("Find α = $(α)")
        if α < 0
            @info "Wrong curvature, α is negative"
            α_t *= 2.0
            println("Should continue to next iter")
            continue
        end
        # actual step
        do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)
        E_t2 = do_compute_energy(Ham, psiks)
        println("Trial energy 2: E_t2 = $(E_t2)")
        if E_t2 > E_old
            @warn "Energy is not decreasing"
            α_t *= 0.1
        else
            println("Accept α=$(α)")
            IS_SUCCESS = true
            break
        end
    end

    if IS_SUCCESS
        return α
    else
        return α_safe
    end

end