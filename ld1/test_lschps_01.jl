using Printf

include("RadialGrid.jl")
include("starting_potential.jl")
include("lschps.jl")

function main()

    Zval = 14.0
    Zed = Zval # no need to convert it to Ry, how about the sign?
    Nspin = 1
    Nwf = 5
    nn = [1, 2, 2, 3, 3]
    ll = [0, 0, 1, 0, 1] 
    oc = [2.0, 2.0, 6.0, 2.0, 2.0]
    enl = zeros(Float64, Nwf)

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf

    rmax = 100.0
    xmin = -8.0 # iswitch = 1
    dx = 0.008 # iswitch = 1
    ibound = false # default

    # Initialize radial grid
    grid = RadialGrid(rmax, Zval, xmin, dx, ibound)
    @info "Nrmesh = $(grid.Nrmesh)"

    Nrmesh = grid.Nrmesh
    v0 = zeros(Float64, Nrmesh)
    vxt = zeros(Float64, Nrmesh)
    vpot = zeros(Float64, Nrmesh, 2)
    enne = 0.0
    starting_potential!(
        Nrmesh, Zval, Zed,
        Nwf, oc, nn, ll,
        grid.r, enl, v0, vxt, vpot, enne, Nspin
    )

    println("After starting_potential:")
    println("v0   = ", v0[1:2])
    println("vxt  = ", vxt[1:2])
    println("vpot1 = ", vpot[1:2,1])
    println("vpot2 = ", vpot[1:2,2])
    println("enl = ", enl[1:Nwf])

    # Solve for all states
    thresh0 = 1.0e-10
    psi = zeros(Float64, Nrmesh, Nwf)
    nstop = 0
    mode = 1
    for iwf in 1:Nwf
        @views psi1 = psi[:,iwf] # zeros wavefunction
        enl[iwf], nstop = lschps!( mode, Zval, thresh0, 
            grid, nn[iwf], ll[iwf], enl[iwf], vpot, psi1
        )
    end

    for iwf in 1:Nwf
        println("outside lschps: enl = ", enl[iwf])
        # println("psi[1] = ", psi[1,iwf])
    end
    println("Pass here")


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
