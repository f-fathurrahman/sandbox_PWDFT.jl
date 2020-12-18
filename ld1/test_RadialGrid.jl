using Printf
include("RadialGrid.jl")

function main()
    rmax = 100.0
    zval = 14.0
    xmin = -7.0 # iswitch = 1
    dx = 0.008
    ibound = false

    grid = RadialGrid(10)
    grid = RadialGrid(rmax, zval, xmin, dx, ibound)

    @printf("Nrmesh = %d\n", grid.Nrmesh)
    @printf("r = %18.10f\n", grid.r[10])

    println("Pass here ...")
end

main()