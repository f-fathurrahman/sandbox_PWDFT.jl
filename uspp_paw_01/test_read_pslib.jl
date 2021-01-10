using Printf
import LightXML

include("PsPot_UPF.jl")

function test_read_pslib()
    #psp = PsPot_UPF("/home/efefer/pseudo/PSLIB/Pd.pbe-n-rrkjus_psl.0.3.0.UPF")
    #psp = PsPot_UPF("/home/efefer/pseudo/PSLIB/H.pbe-rrkjus_psl.0.1.UPF")
    psp = PsPot_UPF("/home/efefer/pseudo/PSLIB/H.pbe-kjpaw_psl.0.1.UPF")
    println("Nproj = ", psp.Nproj)
    println("is_nlcc = ", psp.is_nlcc)
    println("is_ultrasoft = ", psp.is_ultrasoft)
    println("is_paw = ", psp.is_paw)
end

test_read_pslib()
