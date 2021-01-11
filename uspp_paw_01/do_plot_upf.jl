using Printf
import LightXML

import PyPlot
const plt = PyPlot
plt.matplotlib.style.use("dark_background")

include("PsPot_UPF.jl")

function do_plot(filename::String)
    println("filename = ", filename)
    psp = PsPot_UPF(filename)

    #plt.clf()
    #plt.plot(psp.r, psp.V_local)
    #plt.xlim(0.0, 5.0)
    #plt.title("V local of " * filename)
    #plt.savefig("IMG_Vlocal.pdf")

    Nq = size(psp.qfuncl,2)
    nqlc = size(psp.qfuncl,3)
    println("Nproj = ", psp.Nproj)
    println("proj_l = ", psp.proj_l)
    println("Nq = ", Nq)
    println("nqlc = ", nqlc)

    Nproj = psp.Nproj
    proj_l = psp.proj_l
    for nb in 1:Nproj
        plt.clf()
        for mb in nb:Nproj
            ijv = round(Int64, mb*(mb-1)/2) + nb
            l1 = proj_l[nb]
            l2 = proj_l[mb]
            for l in range( abs(l1-l2), stop=(l1+l2), step=2)
                label_txt = string(mb)*"_"*string(nb)*"_"*string(l)
                plt.plot(psp.r, psp.qfuncl[:,ijv,l+1], label=label_txt)
            end
        end
        plt.xlim(0.0, 5.0)
        plt.title("qfuncl of "*psp.atsymb)
        plt.legend()
        filesave = "IMG_qfuncl_nb_"*string(nb)*".png"
        plt.savefig(filesave, dpi=150)
    end

end

@assert length(ARGS) > 0
do_plot(ARGS[1])

