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
