import Serialization

import PyPlot
const plt = PyPlot

function main()

    wgrid = Serialization.deserialize("wgrid.data")    
    εr = Serialization.deserialize("epsr.data")
    εi = Serialization.deserialize("epsi.data")

    println("εr[1,1] = ", εr[1,1])

    plt.clf()
    plt.plot(wgrid, εr[1,:], label="x-real")
    plt.plot(wgrid, εi[1,:], label="x-imag")
    plt.grid(true)
    plt.xlim(0,15)
    plt.savefig("IMG_epsilon.pdf")
end

main()