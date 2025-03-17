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

# %%
Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");
# Compute this once and for all
Ham.energies.NN = calc_E_NN(Ham.atoms);

# %% [markdown]
# We need to set some parameters manually:

# %%
use_smearing = false
kT = 0.0
if pwinput.occupations == "smearing"
    use_smearing = true
    kT = pwinput.degauss*0.5 # convert from Ry to Ha
    Ham.electrons.kT = kT
end
# No need to set kT here, it is already default to 0

# XXX This is probably not needed. We are not directly initializing Rhoe.
# XXX This might be useful for SCF which need this parameter to be
# XXX set for spin-polarized calculations 
if pwinput.nspin == 2
    starting_magnetization = pwinput.starting_magnetization
else
    starting_magnetization = nothing
end

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
Random.seed!(1234)
# This will take into account whether the overlap operator is needed or not
psiks = rand_BlochWavefunc(Ham);

# Prepare Haux (random numbers)

#=
# For Haux, choose between generic symmetric Haux:
Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
for ikspin in 1:Nkspin
    Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
    Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' )
end
=#

# or diagonal Haux:
# Prepare Haux (random numbers)
Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
for ikspin in 1:Nkspin
    Haux[ikspin] = diagm(0 => sort(randn(Float64, Nstates)))
end

# %%
Haux[1]

# %% [markdown]
#
# Initialize gradients and related quantities:

# %%
# Gradients, subspace Hamiltonian
g = zeros_BlochWavefunc(Ham)
Kg = zeros_BlochWavefunc(Ham)
d = zeros_BlochWavefunc(Ham)
#
Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
d_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
for ikspin in 1:Nkspin
    Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    d_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
end

rots_cache = RotationsCache(Nkspin, Nstates);

# %%
Haux_orig = copy(Haux);

# psiks is already orthonormal
# Make Haux diagonal and rotate psiks
# Ham.electrons.ebands are updated here
transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache,
    do_ortho_psi=false, overwrite_Haux=true
)

# %%
Haux[1]

# %%
Haux_orig[1]

# %% [markdown]
# Update Hamiltonian, compute energy and gradients at current psiks and Haux:

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
α = 1.0
do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)
do_step_psiks_Haux!(α, Ham, psiks, Haux, d, d_Haux, rots_cache)

# %%
E1

# %%
do_compute_energy(Ham, psiks)

# %% [markdown]
# Do line minimization:

# %%
E_new, is_success = linmin_quad_v01!(Ham, psiks, Haux, Hsub, g, g_Haux, d, d_Haux, E1)
rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

# %% [markdown]
# Continue iterations

# %%
println("E_new = $(E_new), ΔE = $(E_new - E1)")

# %%
E1 = E_new
