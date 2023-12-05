using Serialization: deserialize

import Plots
const plt = Plots

import PlotThemes
plt.theme(:dark)

r = deserialize("r.dat")
rho_lm = deserialize("rho_lm.dat")

plt.plot(title="Plot", dpi=200)
l2 = size(rho_lm, 2)
for lm in 1:l2
    s = sum(abs.(rho_lm[:,lm,1]))
    println("lm = $lm s = $s")
    if s > 0.1 
        plt.plot!(r, rho_lm[:,lm,1], label="lm="*string(lm))
    end
end
plt.xlims!(0.0, 1.0)
plt.savefig("IMG_rho_lm.png")
