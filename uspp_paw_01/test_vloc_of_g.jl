import LightXML
using Printf
using SpecialFunctions: erf

import PyPlot
const plt = PyPlot

using PWDFT
const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("integ_simpson.jl")
include("PsPot_UPF.jl")
include("vloc_of_g.jl")

function main()

    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """, in_bohr=true, LatVecs = gen_lattice_sc(10.0)
    )

    psp = PsPot_GTH(joinpath(DIR_PSP,"H-q1.gth"))
    psp_upf = PsPot_UPF("/home/efefer/pseudo/HGH/H.pz-hgh.UPF")

    println(psp)
    println(psp_upf)

    pw = PWGrid(15.0, atoms.LatVecs)

    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)

    #Ngl = 101
    #G2_shells = collect(range(0.0, stop=10.0, length=Ngl))

    Vgl = zeros(Float64, Ngl)
    for igl in 1:Ngl
        Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
    end

    Vgl_upf = init_Vloc_G(
        psp_upf.r, psp_upf.rab,
        psp_upf.V_local, psp_upf.zval,
        Ngl, G2_shells,
        pw.CellVolume
    )

    for igl in 1:Ngl
        @printf("%3d %18.10f %18.10f\n", igl, Vgl[igl], Vgl_upf[igl])
    end


    strf = calc_strfact( atoms, pw )
    Npoints = prod(pw.Ns)
    Ng = pw.gvec.Ng
    idx_g2shells = pw.gvec.idx_g2shells
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume

    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)

    V_Ps_loc_upf = zeros(Float64, Npoints)
    Vg_upf = zeros(ComplexF64, Npoints)

    isp = 1
    for ig in 1:Ng
        ip = idx_g2r[ig]
        igl = idx_g2shells[ig]
        Vg[ip] = strf[ig,isp] * Vgl[igl]
        Vg_upf[ip] = strf[ig,isp] * Vgl_upf[igl]
    end
    #
    V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints / CellVolume
    #
    V_Ps_loc_upf[:] = V_Ps_loc_upf[:] + real( G_to_R(pw, Vg_upf) ) * Npoints / CellVolume

    println("avg V_Ps_loc     = ", sum(V_Ps_loc)/Npoints)
    println("avg V_Ps_loc_upf = ", sum(V_Ps_loc_upf)/Npoints)

    for ip in 1:50
        @printf("%3d %18.10f %18.10f %18.10e\n", ip,
            V_Ps_loc[ip], V_Ps_loc_upf[ip], abs(V_Ps_loc[ip]-V_Ps_loc_upf[ip]))
    end

    plt.clf()
    plt.plot(sqrt.(G2_shells), Vgl, marker="o", label="Vgl")
    plt.plot(sqrt.(G2_shells), Vgl_upf, marker="o", label="Vgl_upf")
    plt.legend()
    plt.grid(true)
    plt.xlim(0.0,sqrt(10.0))
    plt.savefig("IMG_Vgl.pdf")
end

main()