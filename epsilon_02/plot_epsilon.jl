import Serialization

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)
plt.rc("font", size=14)

function main()

    wgrid = Serialization.deserialize("wgrid.data")    
    εr = Serialization.deserialize("epsr.data")
    εi = Serialization.deserialize("epsi.data")

    println("εr[1,1] = ", εr[1,1])

    plt.figure(figsize=(7,8))
    plt.clf()
    plt.plot(wgrid, εr[1,:], label="\$\\epsilon_{1}\$")
    plt.plot(wgrid, εi[1,:], label="\$\\epsilon_{2}\$")
    plt.grid(true)
    plt.xlim(0,16)
    plt.xlabel("\$\\hbar\\omega\$ (eV)")
    plt.ylabel("\$\\epsilon\$")
    plt.legend()
    plt.savefig("IMG_epsilon.pdf")
end

main()