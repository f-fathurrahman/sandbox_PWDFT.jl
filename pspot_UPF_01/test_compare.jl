import LightXML
using Printf

using PWDFT
const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("PsPot_UPF.jl")

function main()
    #pspfile1 = "C-q4.gth"
    #pspfile2 = "/home/efefer/pseudo/HGH/C.pz-hgh.UPF"

    #pspfile1 = "Si-q4.gth"
    #pspfile2 = "/home/efefer/pseudo/HGH/Si.pz-hgh.UPF"

    pspfile1 = "Pt-q10.gth"
    pspfile2 = "/home/efefer/pseudo/HGH/Pt.pz-hgh.UPF"

    psp = PsPot_GTH(joinpath(DIR_PSP,pspfile1))
    println(psp)

    psp_upf = PsPot_UPF(pspfile2)
    println(psp_upf)

    idx_r = 10
    r = psp_upf.r[idx_r]
    @printf("r = %18.10f\n", r)
    @printf("V_local = %18.10f\n", psp_upf.V_local[idx_r])
    @printf("eval_V_local = %18.10f\n", PWDFT.eval_Vloc_R(psp, r))

    proj_func = psp_upf.proj_func
    Nproj = psp_upf.Nproj

    channel_idx = zeros(Int64,Nproj) # channel in l, needed for GTH expr
    iprj = 0
    for l in 0:psp.lmax
        for i in 1:psp.Nproj_l[l+1]
            iprj = iprj + 1
            channel_idx[iprj] = i
        end
    end

    println("channel_idx = ", channel_idx)
    println("Nproj_l = ", psp.Nproj_l)

    for iprj in 1:Nproj
        l = psp_upf.proj_l[iprj]
        println("\nl = ", l)

        upf_proj = proj_func[idx_r,iprj]
        println("upf_proj = ", upf_proj)
        
        prj_analytic = PWDFT.eval_proj_R(psp, l, channel_idx[iprj], r)*r
        println("r*eval_proj_R = ", prj_analytic)

        println("ratio = ", prj_analytic/upf_proj)
    end
    println("Pass here")
end

main()