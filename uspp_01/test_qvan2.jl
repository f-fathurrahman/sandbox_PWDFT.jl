using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../ylm_real/Ylm_real_qe.jl")
include("calc_clebsch_gordan.jl")
include("calc_qradG.jl")
include("PsPotNL_UPF.jl")
include("qvan2.jl")

include("../pwscf_02/PWSCFInput.jl")

println("ARGS = ", ARGS)
@assert length(ARGS) == 1
pwinput = PWSCFInput(ARGS[1])

#
# Main program here
#

atoms = pwinput.atoms
#println(atoms)

dual = pwinput.ecutrho/pwinput.ecutwfc
pw = PWGrid(pwinput.ecutwfc, atoms.LatVecs, dual=dual)
#println(pw)

Nspecies = atoms.Nspecies    
pspfiles = pwinput.pspfiles
pspots = Vector{PsPot_UPF}(undef,Nspecies)
for isp in 1:Nspecies
    pspots[isp] = PsPot_UPF( pspfiles[isp] )
    PWDFT._build_prj_interp_table!( pspots[isp], pw )
end

pspotNL = PsPotNL_UPF(atoms, pw, pspots)

Ng = pw.gvec.Ng
G2 = pw.gvec.G2
QfuncG = zeros(ComplexF64,Ng)
println("Ng = ", Ng)

nh = pspotNL.nh
lmaxkb = pspotNL.lmaxkb

# Input
isp = 2
ih = 1
jh = 1

@assert ih <= nh[isp]
@assert jh <= nh[isp]
@assert isp <= Nspecies

# Iinput
lmaxq = 2*lmaxkb + 1 # using 1-indexing
ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
_lmax = lmaxq - 1 # or 2*lmaxkb
Ylm_real_qe!(_lmax, pw.gvec.G, ylmk0) # Ylm_real_qe accept l value starting from 0

qvan2!(pspotNL, ih, jh, isp, G2, ylmk0, QfuncG)
for i in 1:5
    @printf("%5d [%18.10f,%18.10f]\n", i, real(QfuncG[i]), imag(QfuncG[i]))
end

println("sum QfuncG = ", sum(abs.(QfuncG)))