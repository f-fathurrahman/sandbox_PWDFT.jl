import LightXML
import PWDFT: PsPot_UPF

function test_main()
    prefix_dir = "/home/efefer/pseudo/PSLIB/"
    files = readlines("PSLIB_PAW.dat")
    for f in files
        filepath = joinpath(prefix_dir, f)
        println("\nfilepath = ", filepath)
        #
        # Also check PsPot_UPF
        psp = PsPot_UPF(joinpath(prefix_dir, f))
        ##
        println("psp.lmax = ", psp.lmax)
        println("psp.lmax_rho = ", psp.lmax_rho)
        println("psp.nqf = ", psp.nqf)
        println("psp.nqlc = ", psp.nqlc)
        println("psp.q_with_l = ", psp.q_with_l)
        println("psp.is_paw = ", psp.is_paw)
        println("psp.dx = ", psp.dx)
        println("psp.zmesh = ", psp.zmesh)
    end
end

test_main()
