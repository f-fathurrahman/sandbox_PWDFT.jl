
apwfr = apwlo_vars.apwfr;
lofr = apwlo_vars.lofr;

# COMPLEX(8) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)

#= apwfr(1:nrmtmax, 1:2, 1:apwordmax, 0:lmaxapw, 1:natmtot)
Plots.plot(apwfr_elk[:,1,1,0,1]); Plots.plot!(apwfr[1][0][1][:,1])
Plots.plot(apwfr_elk[:,1,1,1,1]); Plots.plot!(apwfr[1][1][1][:,1])
Plots.plot(apwfr_elk[:,1,1,2,1]); Plots.plot!(apwfr[1][2][1][:,1])
Plots.plot(apwfr_elk[:,1,1,3,1]); Plots.plot!(apwfr[1][3][1][:,1])
Plots.plot(apwfr_elk[:,1,1,4,1]); Plots.plot!(apwfr[1][4][1][:,1])
Plots.plot(apwfr_elk[:,1,1,5,1]); Plots.plot!(apwfr[1][5][1][:,1])
Plots.plot(apwfr_elk[:,1,1,6,1]); Plots.plot!(apwfr[1][6][1][:,1])
Plots.plot(apwfr_elk[:,1,1,7,1]); Plots.plot!(apwfr[1][7][1][:,1])
Plots.plot(apwfr_elk[:,1,1,8,1]); Plots.plot!(apwfr[1][8][1][:,1])
=#

# lofr[1][1]


#=
julia> Plots.plot!(lofr[1][1][:,1])
julia> lofr[1] |> size
(6,)

julia> lofr[1][1] |> size
(497, 2)

Plots.plot(lofr_elk[:,1,1,1]); Plots.plot!(lofr[1][1][:,1])

Plots.plot(lofr_elk[:,1,2,1]); Plots.plot!(lofr[1][2][:,1])

Plots.plot(lofr_elk[:,1,3,1]); Plots.plot!(lofr[1][3][:,1])

Plots.plot(lofr_elk[:,1,4,1]); Plots.plot!(lofr[1][4][:,1])

Plots.plot(lofr_elk[:,1,5,1]); Plots.plot!(lofr[1][5][:,1])

Plots.plot(lofr_elk[:,1,6,1]); Plots.plot!(lofr[1][6][:,1])

size(lofr_elk)
(497, 2, 6, 1)
=#