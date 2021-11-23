using Printf
import LightXML
using PWDFT: PsPot_UPF

include("all_oncv.jl")

function test_ONCV_PBE()
    list_file = split(FILELIST_ONCV_PBE, keepempty=false)
    for f in list_file
        psp = PsPot_UPF(joinpath("/home/efefer/pseudo/ONCV_PBE_v1.2", f))
        println("\nf = ", f)
        println("Nproj = ", psp.Nproj)
        println("proj_l = ", psp.proj_l)
        println("Dion = ")
        display(psp.Dion); println()
        if psp.is_nlcc
            println("Using NLCC")
        end
    end
end

test_ONCV_PBE()
