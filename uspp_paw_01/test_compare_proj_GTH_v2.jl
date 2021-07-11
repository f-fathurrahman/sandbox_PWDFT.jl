using Printf
using SpecialFunctions: sphericalbesselj
import LightXML

using PWDFT
const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("PsPot_UPF.jl")
include("integ_simpson.jl")


function build_interp_proj(psp, ecutwfc, CellVolume, cell_factor)

    dq = 0.01
    ndm = psp.kkbeta
    Nproj = psp.Nproj

    nqx = floor( Int64, (sqrt(2*ecutwfc)/dq + 4)*cell_factor )
    println("nqx = ", nqx)

    interp_table = zeros(Float64,nqx,Nproj)

    aux = zeros(Float64, ndm)
    pref = 4*pi/sqrt(CellVolume)

    for ibeta in 1:Nproj
        l = psp.proj_l[ibeta]
        for iq in 1:nqx
            qi = (iq - 1) * dq
            for ir in 1:psp.kkbeta
                jlqr = sphericalbesselj(l, qi*psp.r[ir])
                aux[ir] = psp.proj_func[ir,ibeta] * psp.r[ir] * jlqr
            end
            vqint = integ_simpson( psp.kkbeta, aux, psp.rab )
            interp_table[iq, ibeta] = vqint * pref
        end
    end

    return interp_table
end


function eval_proj_interp(tab, ibeta, Gm)
    dq = 0.01
    # Interpolation procedure
    px = Gm/dq - floor(Int64, Gm/dq)
    ux = 1.0 - px
    vx = 2.0 - px
    wx = 3.0 - px
    i0 = floor(Int64, Gm/dq) + 1
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    Vq = tab[i0,ibeta] * ux * vx * wx / 6.0 +
         tab[i1,ibeta] * px * vx * wx / 2.0 -
         tab[i2,ibeta] * px * ux * wx / 2.0 +
         tab[i3,ibeta] * px * ux * vx / 6.0
    return Vq
end


function main()
    
    cell_factor = 1.0
    ecutwfc = 15.0
    CellVolume = 16.0^3

    psp_gth = PsPot_GTH(joinpath(DIR_PSP, "Pt-q18.gth"))

    psp_gthnum = PsPot_UPF("/home/efefer/pseudo/HGH/Pt.pz-sp-hgh.UPF")
    interp_table = build_interp_proj(psp_gthnum, ecutwfc, CellVolume, cell_factor)

    ibeta = 1
    Gm = 1.0

    println("Check the ordering of iprjl: ")

    for iprjl in 1:psp_gthnum.Nproj
        Vq = eval_proj_interp(interp_table, iprjl, Gm)
        @printf("%3d %18.10f\n", iprjl, Vq)
    end

    for l in 0:psp_gth.lmax
        for iprj in 1:psp_gth.Nproj_l[l+1]
            Vq = eval_proj_G(psp_gth, l, iprj, Gm, CellVolume)
            @printf("%3d %3d %18.10f\n", l, iprj, Vq)
        end
    end

    iprjl = 0
    for l in 0:psp_gth.lmax
        for iprj in 1:psp_gth.Nproj_l[l+1]
            iprjl = iprjl + 1
            Vq = eval_proj_interp(interp_table, iprjl, Gm)
            @printf("%3d %18.10f\n", iprjl, Vq)
        end
    end

end

main()
