import LightXML
using PWDFT: PsPot_UPF

function test_main()
    filename = ARGS[1]
    println("Filename = ", filename)
    psp = PsPot_UPF(filename)
    println("psp.q_with_l = ", psp.q_with_l)
end

test_main()
