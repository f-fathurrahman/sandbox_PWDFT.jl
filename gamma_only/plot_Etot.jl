using DelimitedFiles

import PyPlot
const plt = PyPlot

function main(prefix::String)
    data = readdlm(prefix)
    Ndata = size(data,1)
    @views Etot = data[2:Ndata,2]
    Etot0 = Etot[1]
    plt.clf()
    plt.plot(data[2:Ndata,1], Etot .- Etot0, marker="o")
    plt.grid()
    plt.title(prefix)
    plt.savefig("IMG_"*prefix*".pdf")
end

@assert length(ARGS) == 1
main(ARGS[1])
