using PWDFT

#include("PsPot_UPF.jl")

function test_main()
    #psp = PsPot_UPF("GBRV_LDA/n_lda_v1.2.uspp.F.UPF2")
    psp = PsPot_UPF("GBRV_LDA/h_lda_v1.4.uspp.F.UPF2")
end

test_main()