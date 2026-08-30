using Printf

#=
import PyPlot
const plt = PyPlot
=#

include("RadialGrid.jl")
include("starting_potential.jl")
include("start_scheq.jl")
include("ascheq.jl")

function test_ascheq_01()

    Zval = 14.0
    Zed = Zval
    Nspin = 1
    Nwf = 5
    nn = [1, 2, 2, 3, 3]
    ll = [0, 0, 1, 0, 1] 
    oc = [2.0, 2.0, 6.0, 2.0, 2.0]

#=
    Zval = 1.0
    Zed = Zval
    Nspin = 1
    Nwf = 1
    nn = [1]
    ll = [0] 
    oc = [1.0]
=#

    enl = zeros(Float64, Nwf)

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf

    rmax = 100.0
    xmin = -7.0 # iswitch = 1
    dx = 0.008 # iswitch = 1
    ibound = false # default

    # Initialize radial grid
    grid = RadialGrid(rmax, Zval, xmin, dx, ibound)

    Nrmesh = grid.Nrmesh
    v0 = zeros(Float64, Nrmesh)
    vxt = zeros(Float64, Nrmesh)
    vpot = zeros(Float64, Nrmesh, 2)
    enne = 0.0
    starting_potential!(
        Nrmesh, Zval, Zed,
        Nwf, oc, nn, ll,
        grid.r, enl, v0, vxt, vpot, noscf = true
    )
    println("After starting_potential:")
    println("v0   = ", v0[1:2])
    println("vxt  = ", vxt[1:2])
    println("vpot1 = ", vpot[1:2,1])
    println("vpot2 = ", vpot[1:2,2])
    println("enl = ", enl[1:Nwf])

    #@. vpot = -Zval/grid.r

    # Solve for all states
    ze2 = -Zval # should be 2*Zval in Ry unit
    thresh0 = 1.0e-10
    psi = zeros(Float64, Nrmesh, Nwf)
    nstop = 0
    iwf = 3
    for iwf in 1:Nwf
        println("\nStart iwf = ", iwf)
        @views psi1 = psi[:,iwf] # zeros wavefunction
        enl[iwf], nstop = ascheq!( nn[iwf], ll[iwf], enl[iwf], grid, vpot, ze2, thresh0, psi1 )
    end

    for iwf in 1:Nwf
        println("outside ascheq: enl = ", enl[iwf])
        # println("psi[1] = ", psi[1,iwf])
    end

#=
    plt.clf()
    for iwf in 1:Nwf
        label_str = "psi-" * string(nn[iwf])  * '-' *string(ll[iwf])
        plt.plot(grid.r, psi[:,iwf], label=label_str)
    end
    plt.xlim(0.0, 3.0)
    plt.grid(true)
    plt.legend()
    plt.savefig("IMG_psi1.png", dpi=150)

    plt.clf()
    plt.plot(grid.r, vpot)
    plt.xlim(0.0, 0.2) # the potential is very localized
    plt.grid(true)
    plt.savefig("IMG_vpot.png", dpi=150)
=#


end

main()
