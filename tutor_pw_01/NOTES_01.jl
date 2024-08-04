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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# Lattice vectors: $\mathbf{a}_{1}$, $\mathbf{a}_{2}$, and $\mathbf{a}_{3}$
# combined into $3\times3$ matrix $\mathbf{h}$ as follows:
# $$
# \mathbf{h}=[\mathbf{a}_{1},\mathbf{a}_{2},\mathbf{a}_{3}]
# $$
# where each lattice vectors are given by columns of $\mathbf{h}$.

# %% [markdown]
# Volume of the cell $\Omega=\det\mathbf{h}$

# %% [markdown]
# Scaled coordinates: $\mathbf{r}=\mathbf{h}\mathbf{s}$
#
# Metric tensor: $\mathcal{G}=\mathbf{h}^{\mathsf{T}}\mathbf{h}$
#
# The distances in scaled coordinates $\mathbf{s}$ related to distances in real coordinates $\mathbf{r}$ by metric tensor as follows
# $$
# \left(\mathbf{r}_{i}-\mathbf{r}_{j}\right)^{2}=\left(\mathbf{s}_{i}-\mathbf{s}_{j}\right)^{\mathsf{T}}\mathcal{G}\left(\mathbf{s}_{i}-\mathbf{s}_{j}\right)
# $$

# %% [markdown]
# Periodic boundary condition can be enforced by
# $$
# \mathbf{r}_{\mathrm{pbc}}=\mathbf{r}-\mathbf{h}\left[\mathbf{h}^{-1}\mathbf{r}\right]_{\mathrm{NINT}}
# $$
# where $\left[\cdots\right]_{\mathrm{NINT}}$ denotes the nearest integer value so that the coordinates $\mathbf{r}_{\mathrm{pbc}}$ will always be within the box centered around the origin of the coordinates system.

# %% [markdown]
# Reciprocal lattice vectors $\mathbf{b}_{i}$ are defined as $\mathbf{b}_{i}\cdot\mathbf{a}_{j}=2\pi\delta_{ij}$. They can be arranged into $3\times3$ matrix:
# $$
# \left[\mathbf{b}_{1},\mathbf{b}_{2},\mathbf{b}_{3}\right]=2\pi\left(\mathbf{h}^{\mathsf{T}}\right)^{-1}
# $$

# %% [markdown]
# Planewave basis expansion

# %% [markdown]
# $$
# f_{\mathbf{G}}^{\mathrm{PW}}(\mathbf{r}) =
# \frac{1}{\sqrt{\Omega}}\exp\left[\imath\mathbf{G}\cdot\mathbf{r}\right] =
# \frac{1}{\sqrt{\Omega}}\exp\left[\imath2\pi\mathbf{g}\cdot\mathbf{s}\right]
# $$

# %% [markdown]
# with reciprocal space vectors
# $$
# \mathbf{G}=2\pi\left(\mathbf{h}^{\mathsf{T}}\right)^{-1}\mathbf{g}
# $$
# where $\mathbf{g}=[j,k,l]$ is a triple of integer values (Miller index).

# %% [markdown]
# Test

# %%
using Printf
using LinearAlgebra
using PWDFT

# %% [markdown]
# $$
# f(G) = \sqrt{8\pi^{3}}r_{\mathrm{loc}}^{3}\exp\left[-\frac{1}{2}G^{2}r_{\mathrm{loc}}^{2}\right]
# $$

# %%
function func_G_analytic(Gl::Float64; rloc=1.5)
    pre2 = sqrt(8*π^3)*rloc^3
    expGr2 = exp(-0.5*(Gl*rloc)^2)
    return pre2*expGr2
end;

# %% [markdown]
# $$
# f(r) = \exp\left[ -\frac{1}{2} \left( \frac{r}{r_{\mathrm{loc}}} \right)^2 \right]
# $$

# %%
function func_R_analytic(rl::Float64; rloc=1.5)
    return exp(-0.5*(rl/rloc)^2)
end;

# %% [markdown]
# Note that the function we choose only depend on the length of G-vector (for fG) and only on the distance between points and the center (we have implicitly define the above center at 0).

# %% [markdown]
# Let's define plane wave basis set in a simple cubic lattice.

# %%
L = 10.0 # bohr
ecutwfc = 5.0 # hartree
# Using higher cutoff will give better results

LatVecs = gen_lattice_sc(L)
pw = PWGrid(ecutwfc, LatVecs)

Ng = pw.gvec.Ng
G2 = pw.gvec.G2
G = pw.gvec.G
CellVolume = pw.CellVolume;

# %% [markdown]
# Let's evaluate $f(G)$, the coefficients of the expansions, using analytic expression.

# %%
fG = zeros(Float64, Ng)
for ig in 1:Ng
    Gl = sqrt(G2[ig])
    fG[ig] = func_G_analytic(Gl)/CellVolume
end;

# %% [markdown]
# Note that we also include `1/CellVolume` factor in the definition of expansion coefficients.
# If we don't include  `1/CellVolume` factor in the definition, we can include it in the definition of expansion of $f(r)$.

# %%
Npoints = prod(pw.Ns)
Rgrid = PWDFT.init_grid_R(pw.Ns, LatVecs);
ip = 10
r = Rgrid[:,ip]
s = 0.0 + im*0.0
for ig in 1:Ng
    Gr = dot(G[:,ig], r)
    s += fG[ig] * exp(im*Gr)
end
fR = real(s)
println("PW expansion = ", fR)

# %%
rl = norm(r)
fR_a = func_R_analytic(rl)
println("analytic = ", fR_a)
println("Difference = ", abs(fR - fR_a))

# %% [markdown]
# Our real space function might have some overlaps with its neighbor, which might explain the difference we see above.

# %%
v1 = LatVecs[:,1]
v2 = LatVecs[:,2]
v3 = LatVecs[:,3]

# Add nearest neighbor contributions
# Result will depend on rloc and CellSize
NumNeighbors = 0
for i in -1:1, j in -1:1, k in -1:1
    if i == 0 && j == 0 && k == 0
        continue
    end
    NumNeighbors += 1
    rn = r + i*v1 + j*v2 + k*v3
    rl = norm(rn)
    fR_a += func_R_analytic(rl)
end
println("NumNeighbors = ", NumNeighbors) # should be 27 - 1 = 26

println("analytic (adding nearest neighbors contrib) = ", fR_a)
println("Difference = ", abs(fR - fR_a))


# %% [markdown]
# # Make into functions

# %%
function test_pw_expansions(pw, func_R, func_G; offset_neighbors=[1,1,1], r=nothing, ip=10)
    LatVecs = pw.LatVecs
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    G = pw.gvec.G
    CellVolume = pw.CellVolume
    println("Ng = ", Ng)
    fG = zeros(Float64, Ng)
    for ig in 1:Ng
        Gl = sqrt(G2[ig])
        fG[ig] = func_G(Gl)/CellVolume # we include 1/CellVolume factor here
    end
    Npoints = prod(pw.Ns)
    Rgrid = PWDFT.init_grid_R(pw.Ns, LatVecs);
    if isnothing(r)
        @assert ip <= Npoints
        r = Rgrid[:,ip]
    end
    println("r = ", r)
    s = 0.0 + im*0.0
    for ig in 1:Ng
        Gr = dot(G[:,ig], r)
        s += fG[ig] * exp(im*Gr)
    end
    fR = real(s)
    println("sum = ", s) # check the complex part
    println("PW expansion = ", fR)
    #
    # real space evalution
    rl = norm(r) # should include
    fR_a = func_R(rl)
    println("Real space (no neighbors) = ", fR_a)
    println("Difference = ", abs(fR - fR_a))
    #
    v1 = LatVecs[:,1]
    v2 = LatVecs[:,2]
    v3 = LatVecs[:,3]
    fR_a_orig = copy(fR_a)
    nnx = offset_neighbors[1]; @assert nnx > 0
    nny = offset_neighbors[2]; @assert nny > 0
    nnz = offset_neighbors[3]; @assert nnz > 0
    num_neighbors = 0
    for i in -nnx:nnx, j in -nny:nny, k in -nnz:nnz
        if i == 0 && j == 0 && k == 0
            continue
        end
        num_neighbors += 1
        rn = r + i*v1 + j*v2 + k*v3
        rl = norm(rn)
        fR_a += func_R(rl)
    end
    println("num_neighbors = ", num_neighbors) # should be 27 - 1 = 26
    println("Analytic (adding nearest neighbors contrib) = ", fR_a)
    println("Effect of neighbor = ", abs(fR_a - fR_a_orig))
    println("Difference = ", abs(fR - fR_a))
    #
    return
end;

# %%
L = 10.0 # bohr
ecutwfc = 2.0 # hartree
LatVecs = gen_lattice_sc(L)
pw = PWGrid(ecutwfc, LatVecs);
func_R(r) = func_R_analytic(r, rloc=1.5) # set rloc
func_G(G) = func_G_analytic(G, rloc=1.5) # set rloc
test_pw_expansions(pw, func_R, func_G; offset_neighbors=[1,1,1], r=[0.4, 0.4, 0.4])

# %%
test_pw_expansions(pw, func_R, func_G; offset_neighbors=[2,2,2], r=[0.4, 0.4, 0.4])

# %%
func_R(10.0), func_R(20.0), func_R(30.0)

# %% [markdown]
# Try using larger rloc:

# %%
L = 10.0 # bohr
ecutwfc = 2.0 # hartree
LatVecs = gen_lattice_sc(L)
pw = PWGrid(ecutwfc, LatVecs);
func_R(r) = func_R_analytic(r, rloc=2.5) # set rloc
func_G(G) = func_G_analytic(G, rloc=2.5) # set rloc
test_pw_expansions(pw, func_R, func_G; offset_neighbors=[1,1,1], r=[0.4, 0.4, 0.4])

# %% jupyter={"outputs_hidden": true}
func_R(10.0), func_R(20.0), func_R(30.0)

# %%
L = 10.0 # bohr
ecutwfc = 5.0 # hartree
LatVecs = gen_lattice_sc(L)
pw = PWGrid(ecutwfc, LatVecs);
func_R(r) = func_R_analytic(r, rloc=2.5) # set rloc
func_G(G) = func_G_analytic(G, rloc=2.5) # set rloc
test_pw_expansions(pw, func_R, func_G; offset_neighbors=[1,1,1], r=[0.4, 0.4, 0.4])
