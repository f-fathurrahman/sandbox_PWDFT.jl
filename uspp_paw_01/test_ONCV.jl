using Printf
import LightXML

include("PsPot_UPF.jl")
include("all_oncv.jl")

function test_ONCV_PBE()
    list_file = split(FILELIST_ONCV_PBE, keepempty=false)
    for f in list_file
        psp = PsPot_UPF(joinpath("/home/efefer/pseudo/ONCV_PBE", f))
        println("\nf = ", f)
        println("Nproj = ", psp.Nproj)
        if psp.is_nlcc
            println("Using NLCC")
        end
    end
end

test_ONCV_PBE()
