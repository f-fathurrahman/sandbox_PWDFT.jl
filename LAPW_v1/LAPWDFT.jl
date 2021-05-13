module LAPWDFT

using Printf
using LinearAlgebra
using OffsetArrays

using PWDFT

include("r3frac.jl")
include("r3mv.jl")

export init_zero!

include("LatticeVars.jl")
export LatticeVars

include("AtomicVars.jl")
export AtomicVars

include("AtomicSpeciesVars.jl")
export AtomicSpeciesVars, init_nuclear_pot!

include("MuffinTins.jl")
export MuffinTins, init_packed_mtr!

include("APWLOVars.jl")
export APWLOVars

include("readspecies.jl")
export readspecies!

include("mtdmin.jl")
export mtdmin

include("checkmt.jl")
export checkmt!

include("genrmesh.jl")
export genrmesh!

include("radnucl.jl")
export radnucl

include("potnucl.jl")
export potnucl!

include("XC_funcs/XC_x_slater.jl")
include("XC_funcs/XC_c_vwn.jl")

include("polynm.jl")
export polynm

include("splint.jl")
export splint

include("splintwp.jl")
export splintwp!

include("wsplint.jl")
export wsplint!

include("wsplintp.jl")
export wsplintp

include("solve_atom.jl")
export solve_atom!

include("rdirac.jl")
export rdirac!

include("rdiracint.jl")
export rdiracint!

include("poly4i.jl")
export poly4i

include("poly3.jl")
export poly3

end
