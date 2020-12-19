using Printf

include("RadialGrid.jl")
include("starting_potential.jl")

function main()

    Zval = 14.0
    Zed = Zval
    Nspin = 1
    Nwf = 5
    nn = [1, 2, 2, 3, 3]
    ll = [0, 0, 1, 0, 1]
    oc = [2.0, 2.0, 6.0, 2.0, 2.0]
    enl = zeros(Float64,Nwf)

    rmax = 100.0
    xmin = -7.0 # iswitch = 1
    dx = 0.008
    ibound = false

    grid = RadialGrid(rmax, Zval, xmin, dx, ibound)

    Nrmesh = grid.Nrmesh
    v0 = zeros(Float64, Nrmesh)
    vxt = zeros(Float64, Nrmesh)
    vpot = zeros(Float64, Nrmesh, 2)
    enne = 0.0
    starting_potential!( Nrmesh, Zval, Zed,
        Nwf, oc, nn, ll,
        grid.r, enl, v0, vxt, vpot, enne, Nspin
    )

    println("After starting_potential:")
    println("v0   = ", v0[1:2])
    println("vxt  = ", vxt[1:2])
    println("vpot1 = ", vpot[1:2,1])
    println("vpot2 = ", vpot[1:2,2])
    println("enl = ", enl[1:Nwf])

    println("Pass here ...")
end

main()
