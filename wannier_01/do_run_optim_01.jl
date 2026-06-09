include("../../run_optim.jl")

loc = omega_loc(p,A)
# writedlm("data",loc)
println(maximum(loc))
println("Max loc = ", maximum(loc))

#include("plot_free_wannier.jl")
