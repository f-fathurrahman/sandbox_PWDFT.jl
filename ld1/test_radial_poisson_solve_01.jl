using Revise, Infiltrator
using Printf
import LinearAlgebra

#=
# Use import to avoid potential name clash
import Plots, PlotThemes
Plots.theme(:dark)
=#

using PWDFT


includet("RadialGrid.jl")
includet("starting_potential.jl")
includet("start_scheq.jl")
includet("ascheq.jl")
includet("lschps.jl")
includet("radial_poisson_solve.jl")

function test_radial_poisson_solve()

    Zval = 14.0
    Zed = Zval # no need to convert it to Ry, how about the sign?
    Nspin = 1
    Nwf = 5
    nn = [1, 2, 2, 3, 3]
    ll = [0, 0, 1, 0, 1] 
    oc = [2.0, 2.0, 6.0, 2.0, 2.0]
    enl = zeros(Float64, Nwf)
    # FIXME: define isw: spin index

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf

    rmax = 100.0
    #xmin = -8.0 # iswitch=1, rel=1
    xmin = -7.0 # iswitch=1, rel=0
    dx = 0.008 # iswitch=1, rel=1
    ibound = false # default

    # Initialize radial grid
    grid = RadialGrid(rmax, Zval, xmin, dx, ibound)
    @info "Nrmesh = $(grid.Nrmesh)"

    Nrmesh = grid.Nrmesh
    v0 = zeros(Float64, Nrmesh)
    vxt = zeros(Float64, Nrmesh)
    vpot = zeros(Float64, Nrmesh, 2)
    starting_potential!(
        Nrmesh, Zval, Zed,
        Nwf, oc, nn, ll,
        grid.r, enl, v0, vxt, vpot
    )


    # Solve for all states
    thresh0 = 1.0e-10
    psi = zeros(Float64, Nrmesh, Nwf)
    nstop = 0
    mode = 1 # for lschps
    ze2 = -Zval # should be 2*Zval in Ry unit
    for iwf in 1:Nwf
        @views psi1 = psi[:,iwf] # zeros wavefunction
        #enl[iwf], nstop = lschps!( mode, Zval, thresh0, 
        #    grid, nn[iwf], ll[iwf], enl[iwf], vpot, psi1
        #)
        enl[iwf], nstop = ascheq!(
            nn[iwf], ll[iwf], enl[iwf], grid, vpot, ze2, thresh0, psi1, nstop
        )
    end

    println("Energy levels:")
    for iwf in 1:Nwf
        @printf("%3d %18.10f\n", iwf, enl[iwf])
    end

    #
    # calculate charge density (spherical approximation)
    #
    rho = zeros(Float64, Nrmesh, Nspin)
    ispin = 1
    for iwf in 1:Nwf, ir in 1:Nrmesh
        # this is for ispin=1
        rho[ir,ispin] += oc[iwf] * psi[ir,iwf]^2
    end
    integRho = PWDFT.integ_simpson(Nrmesh, rho[:,1], grid.rab) 
    println("integRho = ", integRho)

    V_h = zeros(Float64, Nrmesh )
    radial_poisson_solve!(0, 2, grid, rho, V_h)

    @infiltrate

    return

end

