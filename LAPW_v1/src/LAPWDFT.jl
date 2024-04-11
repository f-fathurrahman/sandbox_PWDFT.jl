module LAPWDFT

using Printf
using LinearAlgebra
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT: Atoms, PWGrid, KPoints, LibxcXCCalculator,
             SymmetryInfo, R_to_G!, G_to_R!,
             calc_epsxc_Vxc_LDA, calc_epsxc_Vxc_LDA!,
             calc_Vxc_LDA, calc_Vxc_LDA!,
             calc_epsxc_LDA, calc_epsxc_LDA!

using Infiltrator

include("r3frac.jl")
include("r3mv.jl")

export init_zero!

include("LatticeVars.jl")
export LatticeVars

include("AtomicVars.jl")
export AtomicVars

include("AtomicSpeciesVars.jl")
export AtomicSpeciesVars, init_nuclear_pot!

include("SphericalHarmonicTransform.jl")
export SphericalHarmonicTransform

include("MuffinTins.jl")
export MuffinTins, init_packed_mtr!

include("APWLOVars.jl")
export APWLOVars

include("readspecies.jl")
export readspecies!

include("read_elk_input.jl")
export read_elk_input
export create_atoms_from_elk_input

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


# FIXME: Move to PWDFT?
include("Libxc_old.jl")
#include("XCCalculator.jl")
include("LDA_PW92.jl")
#export XCCalculator, LibxcXCCalculator
export calc_epsxc_PW92
export calc_Vxc_PW92, calc_Vxc_PW92!
export calc_epsxc_Vxc_PW92, calc_epsxc_Vxc_PW92!

#include("XC_funcs/XC_x_slater.jl")
#include("XC_funcs/XC_c_vwn.jl")
#include("XC_funcs/XC_c_pw.jl")

include("polynm.jl")
export polynm

include("factnm.jl")
export factnm

include("spline.jl")
export spline!

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

include("allatoms.jl")
export allatoms!

include("calc_sfacg.jl")
export calc_sfacg

include("rhoinit.jl")
export rhoinit!

include("sbessel.jl")
export sbessel!

include("genjlgprmt.jl")
export genjlgprmt!

include("z_to_rf_mt.jl")
include("z_to_rf_lm.jl")
export z_to_rf_lm!, z_to_rf_mt!

include("r_to_zf_mt.jl")
include("r_to_zf_lm.jl")
export r_to_zf_lm!, r_to_zf_mt!

include("rf_mt_c_to_f.jl")
export rf_mt_c_to_f!

include("rf_interp.jl")
export rf_interp!

include("genylmv.jl")
include("genylmg.jl")
include("genrlmv.jl")
export genylmg!, genylmv!, genrlmv!

include("sphcover.jl")
include("sctovec.jl")
export sphcover!, sctovec!

include("zpotclmt.jl")
export zpotclmt!

include("zpotcoul.jl")
export zpotcoul!

include("potcoul.jl")
export potcoul!

include("debug_main.jl")

end

