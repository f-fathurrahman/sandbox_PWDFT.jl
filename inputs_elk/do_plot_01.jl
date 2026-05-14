using Plots, PlotThemes
theme(:dark)

apwfr = apwlo_vars.apwfr;
# for isp=1
rgrid = mt_vars.rlmt[1][:,1]; # r^l where l=1

ia = 1;
#io = 1; # derivative order: io=1, 0-th deriv
l = 1; # from 0:lmaxapw, usually lmaxapw=9
fig = plot(rgrid, apwfr[ia][l][1][:,1], label="ord=1")
plot!(fig, rgrid, apwfr[ia][l][2][:,1], label="ord=2")
