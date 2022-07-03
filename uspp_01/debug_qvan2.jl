using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT

include("../ylm_real/Ylm_real_qe.jl")
include("calc_clebsch_gordan.jl")
include("calc_qradG.jl")
include("PsPotNL_UPF.jl")

include("../pwscf_02/PWSCFInput.jl")

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

#
# Debug qvan2
#

isp = 1
indv = pspotNL.indv
nhtolm = pspotNL.nhtolm
nh = pspotNL.nh
lpx = pspotNL.lpx
lpl = pspotNL.lpl
ap = pspotNL.ap
lmaxkb = pspotNL.lmaxkb
qradG = pspotNL.qradG

Ng = pw.gvec.Ng
G2 = pw.gvec.G2
QfuncG = zeros(ComplexF64,Ng)

# Input
ih = 2
jh = 8
@assert ih <= nh[isp]
@assert jh <= nh[isp]

# Iinput
lmaxq = 2*lmaxkb + 1 # using 1-indexing
ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
_lmax = lmaxq - 1 # or 2*lmaxkb
Ylm_real_qe!(_lmax, pw.gvec.G, ylmk0) # Ylm_real_qe accept l value starting from 0

println("-----------")
println("debug qvan2")
println("-----------")

nb = indv[ih,isp]
mb = indv[jh,isp]

if nb >= mb
    ijv = round(Int64, nb * (nb - 1) / 2 + mb)
else
    ijv = round(Int64, mb * (mb - 1) / 2 + nb)
end

ivl = nhtolm[ih,isp]
jvl = nhtolm[jh,isp]

dq = 0.01 # XXX HARDCODED

@printf("ih,jh = (%3d,%3d) mb,nb = (%3d,%3d) ivl,jvl = (%3d,%3d) ijv = %3d\n", 
    ih, jh, mb, nb, ivl, jvl, ijv)

for lm in 1:lpx[ivl,jvl]
    #
    lp = lpl[ivl,jvl,lm] # combined (l+1) index (using default 1-index based array)

    # Hardcoded assertion
    @assert lp >= 1
    @assert lp <= 49

    # finds angular momentum l corresponding to combined index lp (l is 
    # actually l+1 because this is the way qrad is stored, check init_us_1)
    if lp == 1
        l = 1
        sig = 1.0
        ind = 1 # real
    elseif lp <= 4
        l = 2
        sig = -1.0
        ind = 2 # imag
    elseif lp <= 9
        l = 3
        sig = -1.0
        ind = 1
    elseif lp <= 16
        l = 4
        sig = 1.0
        ind = 2
    elseif lp <= 25
        l = 5
        sig = 1.0
        ind = 1
    elseif lp <= 36
        l = 6
        sig = -1.0
        ind = 2
    else
        l = 7
        sig = -1.0
        ind = 1
    end

    println("lm = ", lm, " lp = ", lp, " l = ", l, " (-im)^l = ", (-im)^(l-1), " sig = ", sig, " ind = ", ind)
    # Physics l

    # sig = sig * ap(lp, ivl, jvl) # access Clebsch-Gordan coef
    prefact = (-im)^l * ap[lp,ivl,jvl]

    for ig in 1:Ng
        # calculates quantites depending on the module of G only when needed
        Gm = sqrt(G2[ig])
        px = Gm/dq - floor(Int64, Gm/dq)
        ux = 1.0 - px
        vx = 2.0 - px
        wx = 3.0 - px
        i0 = floor(Int64, Gm/dq) + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        uvx = ux * vx * (1/6)
        pwx = px * wx * 0.5
        work = qradG[isp][i0,ijv,l] * uvx * wx +
               qradG[isp][i1,ijv,l] * pwx * vx -
               qradG[isp][i2,ijv,l] * pwx * ux +
               qradG[isp][i3,ijv,l] * px * uvx
        QfuncG[ig] = QfuncG[ig] + prefact * ylmk0[ig,lp] * work
    end

end

