# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.11.4
#     language: julia
#     name: julia-1.11
# ---

# %%
using Revise

# %%
using LinearAlgebra
using Printf
using Infiltrator
using Random
using PWDFT

# %% [markdown]
# Trying to use `Revise.jl` here:

# %%
includet("smearing.jl")
includet("occupations.jl")
includet("Lfunc.jl")
includet("gradients_psiks_Haux.jl")
includet("utilities_emin_smearing.jl")
includet("prepare_Ham_various.jl")

# %% [markdown]
# Initializing Hamiltonian:

# %%
Ham = prepare_Ham_from_pwinput("PWINPUT");

# %% [markdown]
# Some shortcuts:

# %%
Nspin = Ham.electrons.Nspin
Nkpt = Ham.pw.gvecw.kpoints.Nkpt
Nkspin = Nkpt*Nspin
Nstates = Ham.electrons.Nstates;

# %% [markdown]
# Initialize electronic variables: `psiks` and `Haux`:

# %%
# Symmetric Haux, but not diagonal
function gen_random_symm_Haux(Ham)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates;
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
    end
    return Haux
end

# or diagonal Haux:
# Prepare Haux (random numbers)
function gen_random_diag_Haux(Ham)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates;
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm(0 => sort(randn(Float64, Nstates)))
    end
    return Haux
end

# %%
Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
for ikspin in 1:Nkspin
    Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
end

# %%
# Calculate Hsub
for ispin in 1:Nspin, ik in 1:Nkpt
    Ham.ispin = ispin
    Ham.ik = ik
    ikspin = ik + (ispin-1)*Nkpt
    Hsub[ikspin][:,:] = psiks[ikspin]' * (Ham * psiks[ikspin])
end

# %% [markdown]
# ### Generate initial point: (psiks and Haux)

# %%
Random.seed!(1234)
psiks = rand_BlochWavefunc(Ham)
Haux = gen_random_symm_Haux(Ham);

# %% [markdown]
# ### Do transformation:

# %%
Haux_orig = copy(Haux);
psiks_orig = copy(psiks);
rots_cache = RotationsCache(Nkspin, Nstates);
# psiks is already orthonormal
# Make Haux diagonal and rotate psiks
# Ham.electrons.ebands are updated here
transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache,
    do_ortho_psi=false, overwrite_Haux=true
)

# %% [markdown] jp-MarkdownHeadingCollapsed=true
# Update Hamiltonian, compute energy and gradients at current psiks and Haux.
#
# Haux is stored in ebands and psiks already rotated in case Haux is not diagonal.

# %% [markdown]
# ### Initialize gradients and related quantities:

# %%
# Gradients
g = zeros_BlochWavefunc(Ham)
Kg = zeros_BlochWavefunc(Ham)
d = zeros_BlochWavefunc(Ham)
#
g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
d_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
for ikspin in 1:Nkspin
    g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    d_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
end

# %%
# Update Hamiltonian before evaluating free energy
update_from_ebands!( Ham )
update_from_wavefunc!( Ham, psiks )
E1 = calc_Lfunc( Ham, psiks )
println("E1 = $(E1)")
#
# Calculate gradients
calc_grad_psiks!(Ham, psiks, g, Hsub)
my_Kprec!(Ham, g, Kg)
calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)

# %%
println("Test grad psiks before rotate: $(2*dot(g, psiks))")
println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")

# %%
rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

# %%
println("Test grad psiks after rotate: $(2*dot(g, psiks))")
println("Test grad Haux after rotate: $(dot(Haux, g_Haux))")

# %%
println("Test grad psiks after rotate: $(2*dot(g, psiks_orig))")
println("Test grad Haux orig after rotate: $(dot(Haux_orig, g_Haux))")

# %% [markdown]
# We set the direction (start of minimization):

# %%
# Set direction
for ikspin in 1:Nkspin
    d[ikspin][:,:] = -Kg[ikspin][:,:]
    d_Haux[ikspin][:,:] = -Kg_Haux[ikspin][:,:]
end
constrain_search_dir!(d, psiks)

# Check direction (this is also done in linmin)
gd = 2*real(dot(g,d)) + real(dot(g_Haux, d_Haux))
println("gd = $(gd)")
if gd > 0
    error("Bad step direction")
end

# %% [markdown]
# Trying some steps (this is only for debugging):

# %%
α = -1.0
do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)
do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)

# %%
E1

# %%
do_compute_energy(Ham, psiks)

# %%
Ham.electrons.Focc
