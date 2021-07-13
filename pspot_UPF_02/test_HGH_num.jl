using Printf
import LightXML

include("PsPot_UPF.jl")

function test_HGH_num()
    psp = PsPot_UPF("/home/efefer/pseudo/HGH/I.pz-hgh.UPF")
    println("Nproj = ", psp.Nproj)
    println("is_nlcc = ", psp.is_nlcc)
    println("is_ultrasoft = ", psp.is_ultrasoft)
    println("is_paw = ", psp.is_paw)
end

test_HGH_num()
