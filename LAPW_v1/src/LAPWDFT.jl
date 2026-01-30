module LAPWDFT

using Printf
using LinearAlgebra
using OffsetArrays
using SpecialFunctions: sphericalbesselj
using Serialization: serialize, deserialize

import PWDFT
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

include("SpeciesInfo.jl")
export SpeciesInfo

include("AtomicSpeciesVars.jl")
export AtomicSpeciesVars, init_nuclear_pot!

include("ElectronicChargesStates.jl")
export ElectronicChargesStates

include("SphericalHarmonicTransform.jl")
export SphericalHarmonicTransform

include("backward_SHT.jl")
include("forward_SHT.jl")
export backward_SHT!, forward_SHT!

include("SymmetryVars.jl")
export SymmetryVars

include("findsymlat.jl")
export findsymlat, findsymlat!

include("findsym.jl")
export findsym!

include("findsymcrys.jl")
export findsymcrys!

include("findsymsite.jl")
export findsymsite!

include("MuffinTins.jl")
export MuffinTins, init_packed_mtr!

include("APWLOVars.jl")
export APWLOVars

include("CoreStatesVars.jl")
export CoreStatesVars

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


# Symmetrization
include("symrfir.jl")
export symrfir!


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

# This is maybe required to get similar result with Elk
# In some subroutines, Elk uses full G-vectors (which are not constrained with G2max)
include("GVectorsFull.jl")
export GVectorsFull

include("gensfacgp.jl")
export gensfacgp!

include("genffacgp.jl")
export genffacgp!

include("gencfun.jl")
export gencfun

include("rhoinit.jl"); export rhoinit!
include("maginit.jl"); export maginit!

include("sbessel.jl"); export sbessel!
include("sbesseldm.jl"); export sbesseldm!

include("genjlgprmt.jl")
export genjlgprmt!

include("z_to_rf_mt.jl")
include("z_to_rf_lm.jl")
export z_to_rf_lm!, z_to_rf_mt!

include("r_to_zf_mt.jl")
include("r_to_zf_lm.jl")
export r_to_zf_lm!, r_to_zf_mt!

include("rf_mt_c_to_f.jl"); export rf_mt_c_to_f!
include("rf_mt_f_to_c.jl"); export rf_mt_f_to_c!

include("rf_interp.jl"); export rf_interp!

include("genylmv.jl"); export genylmg!
include("genylmg.jl"); export genylmv!
include("genrlmv.jl"); export genrlmv!

include("sphcover.jl"); export sphcover!
include("sctovec.jl"); export sctovec!

include("roteuler.jl")
export roteuler!

include("ylmroty.jl")
export ylmroty!

include("rlmrot.jl")
export rlmrot!

include("rotrflm.jl")
export rotrflm!

include("rotrfmt.jl")
export rotrfmt!

include("symrfmt.jl")
export symrfmt!

include("zpotclmt.jl")
export zpotclmt!

include("zpotcoul.jl")
export zpotcoul!

include("potcoul.jl"); export potcoul!
include("potxcmt.jl"); export potxcmt!
include("potxcir.jl"); export potxcir!

include("genvsig.jl"); export genvsig!

include("rf_mt_lm.jl"); export rf_mt_lm!

include("genbs.jl"); export genbs!

include("potks.jl"); export potks!

include("gencore.jl"); export gencore!

include("rschrodint.jl"); export rschrodint!

include("findband.jl"); export findband!

include("linengy.jl"); export linengy!

include("genapwfr.jl"); export genapwfr!

include("genlofr.jl"); export genlofr!

include("olprad.jl"); export olprad!

include("hmlrad.jl"); export hmlrad!

include("factr.jl"); export factr
include("wigner3j.jl"); export wigner3j
include("gaunt.jl"); export gaunt
include("gauntyry.jl"); export gauntyry

include("calc_match_coeffs.jl"); export calc_match_coeffs!, calc_match_coeffs

include("zmctmu.jl"); export zmctmu!
include("hmlaa.jl"); export hmlaa!
include("hmlistl.jl"); export hmlistl!
include("hmlalo.jl"); export hmlalo!
include("hmllolo.jl"); export hmllolo!

include("olpaa.jl"); export olpaa!
include("olpistl.jl"); export olpistl!
include("olpalo.jl"); export olpalo!
include("olplolo.jl"); export olplolo!
include("APWLOIntegrals.jl"); export APWLOIntegrals, calc_apwlo_integrals!

include("zfmtinp.jl"); export zfmtinp
include("gen_eigensystem_2nd.jl"); export gen_eigensystem_2nd!
include("gen_eigensystem.jl"); export gen_eigensystem!

include("wavefmt.jl"); export wavefmt!

include("stheta_fd.jl"); export stheta_fd
include("sdelta_fd.jl"); export sdelta_fd
include("occupy.jl"); export occupy!

include("rhomagk.jl"); export rhomagk!
include("rhomagsh.jl"); export rhomagsh!

include("symrvfmt.jl"); export symrvfmt!
include("symrvfir.jl"); export symrvfir!

include("rhocore.jl"); export rhocore!

include("calc_charge.jl"); export calc_charge!, calc_chargemt!
include("rhonorm.jl"); export rhonorm!
include("calc_mag_moment.jl"); export calc_mag_moment!
include("rhomag.jl"); export rhomag!

include("rf_mt_inner_prod.jl")
include("rf_inner_prod.jl")

include("calc_engynn.jl") # XXX should be renamed?
include("calc_Ekin_core.jl")
include("calc_energy_terms.jl")

include("debug_main.jl")
include("debug_scf_01.jl")

end

