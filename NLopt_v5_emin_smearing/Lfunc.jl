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


function reset_rotations!(rots_cache)
    Nkspin = size(rots_cache.rotPrev, 1)
    Nstates = size(rots_cache.rotPrev[1], 1)
    for ikspin in 1:Nkspin
        rots_cache.rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.Urot[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.UrotC[ikspin] = Matrix(1.0*I(Nstates))
    end
    return
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
function update_from_ebands!(Ham, ebands)

    # NOTE: ebands are assumed to be updated outside this function

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    Nelectrons = Ham.electrons.Nelectrons
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    kT = Ham.electrons.kT
    Nspin = Ham.electrons.Nspin
    @assert kT > 1e-10

    E_fermi, mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT,
        Nkpt, wk
    )

    println("mTS = $(mTS), E_fermi=$(E_fermi)")
    # Set some output
    Ham.electrons.E_fermi = E_fermi
    Ham.energies.mTS = mTS
    println("Ham.energies.mTS = ", Ham.energies.mTS)

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

#=
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
=#


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
    for ispin in 1:Nspin, ik in 1:Nkpt
        # Don't forget to set current index for Hamiltonian
        # Need this in case we need to apply op_S
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        #
        ebands[:,ikspin], Urot[ikspin] = eigen(Hermitian(Haux[ikspin]))
        if overwrite_Haux
            Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
        end
        #
        # XXX Check this again, probably do_ortho_psi should always be true
        if do_ortho_psi
            if Ham.need_overlap
                UrotC[ikspin] = inv(sqrt(psiks[ikspin]' * op_S(Ham, psiks[ikspin])))
            else
                UrotC[ikspin] = inv(sqrt(psiks[ikspin]' * psiks[ikspin]))
            end
        end
        UrotC[ikspin] = UrotC[ikspin]*Urot[ikspin] # extra rotation
        psiks[ikspin] = psiks[ikspin]*UrotC[ikspin]
    end
    #
    for ikspin in 1:Nkspin
        # Save previous
        rotPrev[ikspin] = rotPrev[ikspin] * Urot[ikspin]
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
    # Ham.energies.mTS is computed in update_from_ebands!
    return sum(Ham.energies)
end

#=
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
=#


#=
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
=#

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



function linmin_quad_v01!(
    α_t,
    Ham, psiks, Haux, Hsub, g, g_Haux,
    Kg, Kg_Haux,
    d, d_Haux, rots_cache,
    E_old
)

    gd = 2*real(dot(g,d)) + real(dot(g_Haux, d_Haux))
    println("gd = $(gd)")
    α_prev = 0.0
    if gd > 0
        println("ERROR: Bad step direction")
        return E_old, false, α_prev
    end

    NtryMax = 5

    α_safe =  1e-5 # safe step size
    α = α_safe # declare this outside for loop, set to a "safe" value
    α_t_ReduceFactor = 0.1
    α_t_IncreaseFactor = 3.0
    is_success = false
    for itry in 1:NtryMax
        println("--- Begin itry linmin trial step = $(itry) using α_t=$(α_t)")
        #
        do_step_psiks_Haux!(α_t - α_prev, Ham, psiks, Haux, d, d_Haux, rots_cache)
        # make explicit calls to update_* functions
        #
        α_prev = α_t
        E_t = do_compute_energy(Ham, psiks) # this will update ebands, Focc, and Rhoe
        #
        if !isfinite(E_t)
            α_t *= α_t_ReduceFactor
            println("α_t is reduced to=$(α_t)")
            continue # continue
        end
        # prediciton of step size
        c = ( E_t - (E_old + α_t*gd) ) / α_t^2
        α = -gd/(2*c)
        if α < 0
            println("Wrong curvature, α is negative: E_t=$(E_t), E_old=$(E_old)")
            α_t *= α_t_IncreaseFactor
            println("Trial step will become true step. α_t will be set to $(α_t)")
            # calculate gradients
            calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
            calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
            #@infiltrate
            # return trial energy and status
            is_success = true
            return E_t, is_success, α_t
        end
        break
    end
    println("Find α = $(α)")
    
    # actual step
    for itry in 1:NtryMax
        #
        println("--- Begin itry linmin actual step = $(itry) using α=$(α)")
        #
        do_step_psiks_Haux!(α - α_prev, Ham, psiks, Haux, d, d_Haux, rots_cache)
        α_prev = α
        # calculate energy and gradients
        E_t2 = do_compute_energy(Ham, psiks)
        calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
        #
        println("Actual step energy 2: E_t2 = $(E_t2)")
        #
        if !isfinite(E_t2)
            α *= α_t_ReduceFactor
            println("α is reduced to=$(α)")
        end
        # trial energy is higher, reduce α
        if E_t2 > E_old
            α *= α_t_ReduceFactor
            println("Energy is not decreasing, try do decrease α to $(α)")
            continue # continue iteration
        else
            println("Actual step is successful")
            is_success = true
            return E_t2, is_success, α
        end
    end

    # default is unsuccessful try
    return Inf, false, α

end