using Printf
using SpecialFunctions: sphericalbesselj
import LightXML

include("PsPot_UPF.jl")
include("integ_simpson.jl")

function main()
    
    cell_factor = 1.0
    dq = 0.01
    ecutwfc = 15.0 # Ha
    CellVolume = 16.0^3

    nqx = floor( Int64, (sqrt(2*ecutwfc)/dq + 4)*cell_factor )
    println("nqx = ", nqx)

    pspots = Array{PsPot_UPF,1}(undef,2)
    pspots[1] = PsPot_UPF("/home/efefer/pseudo/HGH/C.pbe-hgh.UPF")
    pspots[2] = PsPot_UPF("/home/efefer/pseudo/HGH/O.pbe-hgh.UPF")

    Nspecies = length(pspots)

    Nproj_max = 0
    ndm = 0
    for isp in 1:Nspecies
        psp = pspots[isp]
        kkbeta = psp.kkbeta
        println("kkbeta = ", kkbeta)
        if ndm < kkbeta
            ndm = kkbeta
        end
        if Nproj_max < psp.Nproj
            Nproj_max = psp.Nproj
        end
    end
    println("Nproj_max = ", Nproj_max)
    println("ndm = ", ndm)

    interp_table = zeros(Float64,nqx,Nproj_max,Nspecies)

    aux = zeros(Float64,ndm)
    pref = 4*pi/sqrt(CellVolume)
    println("pref = ", pref)

    for isp in 1:Nspecies
        psp = pspots[isp]
        for ibeta in 1:psp.Nproj
            l = psp.proj_l[ibeta]
            for iq in 1:nqx
                qi = (iq - 1) * dq
                for ir in 1:psp.kkbeta
                    jlqr = sphericalbesselj(l, qi*psp.r[ir])
                    #@printf("%8d %18.10f %18.10f %18.10f %18.10f\n", ir, qi, psp.r[ir],
                    #    jlqr, 2*psp.proj_func[ir,ibeta])
                    aux[ir] = psp.proj_func[ir,ibeta] * psp.r[ir] * jlqr
                end
                #println("sum(aux) = ", 2*sum(aux[1:psp.kkbeta]))
                vqint = integ_simpson( psp.kkbeta, aux, psp.rab )
                interp_table[iq, ibeta, isp] = vqint * pref
            end
        end
    end

    println(pspots[1].r[1:4])

    ibeta = 1
    isp = 2
    for iq in 1:5
        qi = (iq - 1) * dq
        @printf("%18.10f %18.10f\n", qi, interp_table[iq,ibeta,isp]*2)
    end

    Gm = 2.0 # already multiplied by tpiba in QE
    # Interpolation procedure
    px = Gm/dq - floor(Int64, Gm/dq)
    println("px = ", px)
    ux = 1.0 - px
    vx = 2.0 - px
    wx = 3.0 - px
    i0 = floor(Int64, Gm/dq) + 1
    println("i0 = ", i0)
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    Vq = interp_table[i0,ibeta,isp] * ux * vx * wx / 6.0 +
         interp_table[i1,ibeta,isp] * px * vx * wx / 2.0 -
         interp_table[i2,ibeta,isp] * px * ux * wx / 2.0 +
         interp_table[i3,ibeta,isp] * px * ux * vx / 6.0
    @printf("Vq = %18.10f\n", Vq*2)
end

main()
