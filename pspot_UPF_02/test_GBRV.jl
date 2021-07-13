using Printf
import LightXML

include("PsPot_UPF.jl")
include("all_gbrv_files.jl")

function test_GBRV_LDA()
    list_file = split(FILELIST_GBRV_LDA, keepempty=false)
    for f in list_file
        psp = PsPot_UPF(joinpath("./GBRV_LDA", f))
        println()
        println("f = ", f)
        println("Nproj = ", psp.Nproj)
        println("is_nlcc = ", psp.is_nlcc)
        println("is_ultrasoft = ", psp.is_ultrasoft)
        println("is_paw = ", psp.is_paw)
        if !psp.is_nlcc
            @printf("%s does not have NLCC\n", f)
        end
        println("proj_l = ", psp.proj_l)
        println("Dion = ")
        display(psp.Dion); println()
    end
    println("Done for all GBRV_LDA")
end

function test_GBRV_PBE()
    list_file = split(FILELIST_GBRV_PBE, keepempty=false)
    for f in list_file
        psp = PsPot_UPF(joinpath("./GBRV_PBE", f))
        println()
        println("f = ", f)
        println("Nproj = ", psp.Nproj)
        println("is_nlcc = ", psp.is_nlcc)
        println("is_ultrasoft = ", psp.is_ultrasoft)
        println("is_paw = ", psp.is_paw)
        if !psp.is_nlcc
            @printf("%s does not have NLCC\n", f)
        end
        println("proj_l = ", psp.proj_l)
        println("Dion = ")
        display(psp.Dion); println()
    end
    println("Done for all GBRV_PBE")
end

test_GBRV_LDA()
test_GBRV_PBE()
