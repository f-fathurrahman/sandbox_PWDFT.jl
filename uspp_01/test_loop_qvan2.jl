using Printf
using OffsetArrays

using PWDFT

include("create_atoms_N2H4.jl")
include("../ylm_real/Ylm_real_qe.jl")
include("calc_clebsch_gordan.jl")
include("PsPotNL_UPF.jl")

#
# Main program here
#

atoms = create_atoms_N2H4()
#println(atoms)

pw = PWGrid(20.0, atoms.LatVecs, dual=5.0)
#println(pw)

Nspecies = atoms.Nspecies    
pspfiles = [
    "/home/efefer/pseudo/PSLIB/N.pbe-n-rrkjus_psl.0.1.UPF",
    "/home/efefer/pseudo/PSLIB/H.pbe-rrkjus_psl.0.1.UPF"
]

pspots = Vector{PsPot_UPF}(undef,Nspecies)
for isp in 1:Nspecies
    pspots[isp] = PsPot_UPF( pspfiles[isp] )
    PWDFT._build_prj_interp_table!( pspots[isp], pw )
end

pspotNL = PsPotNL_UPF(atoms, pw, pspots)


println("---------------")
println("test loop qvan2")
println("---------------")

isp = 1
indv = pspotNL.indv
nhtolm = pspotNL.nhtolm
nh = pspotNL.nh
lpx = pspotNL.lpx
lpl = pspotNL.lpl

for ih in 1:nh[isp], jh in ih:nh[isp]

    nb = indv[ih,isp]
    mb = indv[jh,isp]

    if nb >= mb
        ijv = round(Int64, nb * (nb - 1) / 2 + mb)
    else
        ijv = round(Int64, mb * (mb - 1) / 2 + nb)
    end

    ivl = nhtolm[ih,isp]
    jvl = nhtolm[jh,isp]

    @printf("ih,jh = (%3d,%3d) mb,nb = (%3d,%3d) ivl,jvl = (%3d,%3d) ijv = %3d\n", 
        ih, jh, mb, nb, ivl, jvl, ijv)

    for lm in 1:lpx[ivl,jvl]
        lp = lpl[ivl,jvl,lm]
        println("lm = ", lm, " lp = ", lp)
        # sig = sig * ap(lp, ivl, jvl) # access Clebsch-Gordan coef
    end
end

