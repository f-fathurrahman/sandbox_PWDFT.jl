include("PsPot_UPF.jl")
include("all_gbrv_files.jl")

function test_GBRV_LDA()
    list_file = split(FILELIST_GBRV_LDA, keepempty=false)
    for f in list_file
        psp = PsPot_UPF(joinpath("./GBRV_LDA", f))
        if !psp.is_nlcc
            @printf("%s does not have NLCC\n", f)
        end
    end
    println("Done for all GBRV_LDA")
end

function test_GBRV_PBE()
    list_file = split(FILELIST_GBRV_PBE, keepempty=false)
    for f in list_file
        psp = PsPot_UPF(joinpath("./GBRV_PBE", f))
        if !psp.is_nlcc
            @printf("%s does not have NLCC\n", f)
        end
    end
    println("Done for all GBRV_PBE")
end

test_GBRV_LDA()
test_GBRV_PBE()
