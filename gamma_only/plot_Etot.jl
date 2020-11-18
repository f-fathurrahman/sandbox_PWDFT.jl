using DelimitedFiles

import PyPlot
const plt = PyPlot

function main(prefix::String)
    Ha2eV = 27.211386018
    data = readdlm(prefix)
    Ndata = size(data,1)
    @views Etot = data[2:Ndata,2]
    Etot0 = Etot[1]
    plt.clf()
    #plt.plot(data[2:Ndata,1], Etot .- Etot0, marker="o")
    plt.plot(data[2:Ndata,1], (Etot .- Etot0)*Ha2eV*1000)
    plt.grid()
    plt.title(prefix)
    plt.ylabel("meV")
    plt.savefig("IMG_"*prefix*".pdf")
end

@assert length(ARGS) == 1
main(ARGS[1])
