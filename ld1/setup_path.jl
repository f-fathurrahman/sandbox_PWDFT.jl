using Revise, Infiltrator
using Printf
import LinearAlgebra

# Use import to avoid potential name clash
import Plots, PlotThemes
Plots.theme(:dark)

using PWDFT

includet("LD1xInput.jl")
includet("RadialGrid.jl")
includet("starting_potential.jl")
includet("start_scheq.jl")
includet("ascheq.jl")
includet("lschps.jl")
includet("radial_poisson_solve.jl")

includet("ld1x_deriv_7pts.jl")
includet("ld1x_deriv2_7pts.jl")
includet("ld1x_find_qi.jl")
