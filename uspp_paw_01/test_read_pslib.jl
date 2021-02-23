using Printf
import LightXML

include("PsPot_UPF.jl")
include("all_pslib_uspp_files.jl")

function test_pslib_uspp()
    list_file = split(FILELIST_PSLIB_USPP, keepempty=false)
    root_pslib = "/home/efefer/pseudo/PSLIB"
    for f in list_file
        psp = PsPot_UPF(joinpath(root_pslib, f))
        println()
        println("f = ", f)
        println("Nproj = ", psp.Nproj)
        println("is_nlcc = ", psp.is_nlcc)
        println("is_ultrasoft = ", psp.is_ultrasoft)
        println("is_paw = ", psp.is_paw)
        if !psp.is_nlcc
            @printf("%s does not have NLCC\n", f)
        end
        println("Dion = ")
        display(psp.Dion); println()
    end
    println("Done for all PSLIB_USPP")
end
test_pslib_uspp()

function test_read_pslib()
    psp = PsPot_UPF("/home/efefer/pseudo/PSLIB/Pd.pbe-n-rrkjus_psl.0.3.0.UPF")
    #psp = PsPot_UPF("/home/efefer/pseudo/PSLIB/H.pbe-rrkjus_psl.0.1.UPF")
    #psp = PsPot_UPF("/home/efefer/pseudo/PSLIB/H.pbe-kjpaw_psl.0.1.UPF")
    println("Nproj = ", psp.Nproj)
    println("is_nlcc = ", psp.is_nlcc)
    println("is_ultrasoft = ", psp.is_ultrasoft)
    println("is_paw = ", psp.is_paw)
end

#test_read_pslib()
