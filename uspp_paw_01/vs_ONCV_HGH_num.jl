using Printf
import LightXML

include("PsPot_UPF.jl")

function main()

    #fnames = [
    #    "/home/efefer/pseudo/ONCV_PBE/I_ONCV_PBE-1.1.upf",
    #    "/home/efefer/pseudo/HGH/I.pbe-hgh.UPF"
    #]

    #fnames = [
    #    "/home/efefer/pseudo/ONCV_PBE/Pt_ONCV_PBE-1.0.upf",
    #    "/home/efefer/pseudo/HGH/Pt.pbe-hgh.UPF"
    #]

    fnames = [
        "/home/efefer/pseudo/ONCV_PBE/Ni_ONCV_PBE-1.0.upf",
        "/home/efefer/pseudo/HGH/Ni.pbe-sp-hgh.UPF"
    ]

    for f in fnames
        psp = PsPot_UPF(f)
        println()
        println(psp.pspfile)
        println("Zval = ", psp.zval)
        println("Nproj = ", psp.Nproj)
        println("proj_l = ", psp.proj_l)
        println("Dion = ")
        display(psp.Dion); println()
        if psp.is_nlcc
            println("Using NLCC")
        else
            println("No NLCC")
        end
    end
end

main()