module LAPWDFT

using Printf
using LinearAlgebra
using OffsetArrays

include("r3frac.jl")
include("r3mv.jl")

include("LatticeVars.jl")
export LatticeVars

include("AtomicVars.jl")
export AtomicVars

include("AtomicSpeciesVars.jl")
export AtomicSpeciesVars

include("MuffinTins.jl")
export MuffinTins

include("APWLOVars.jl")
export APWLOVars

include("readspecies.jl")
export readspecies!

include("mtdmin.jl")
export mtdmin

include("checkmt.jl")
export checkmt

include("genrmesh.jl")
include("radnucl.jl")
include("potnucl.jl")

include("XC_funcs/XC_x_slater.jl")
include("XC_funcs/XC_c_vwn.jl")

include("polynm.jl")
include("splint.jl")
include("splintwp.jl")
include("wsplint.jl")
include("wsplintp.jl")

include("solve_atom.jl")
include("rdirac.jl")
include("rdiracint.jl")
include("poly4i.jl")
include("poly3.jl")

end
