using DelimitedFiles

import PyPlot
const plt = PyPlot

function main()
    plt.clf()

    data = readdlm("ETOT_CO2_extrap2nd.dat")
    Ndata = size(data,1)
    @views Etot = data[2:Ndata,2]
    Etot0 = Etot[1]    
    plt.plot(data[2:Ndata,1], Etot .- Etot0, label="extrap2nd")

    data = readdlm("ETOT_CO2_noextrap.dat")
    Ndata = size(data,1)
    @views Etot = data[2:Ndata,2]
    #Etot0 = Etot[1]    
    plt.plot(data[2:Ndata,1], Etot .- Etot0, label="noextrap")

    plt.legend()
    plt.grid()
    plt.title("Compare")
    plt.savefig("IMG_compare.pdf")
end

main()
