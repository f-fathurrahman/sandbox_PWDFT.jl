using Printf
include("RadialGrid.jl")

function main()

    zval = 14.0 # Depend on the atom

    rmax = 100.0
    xmin = -7.0 # iswitch = 1
    dx = 0.008
    ibound = false

    # grid = RadialGrid(10)
    grid = RadialGrid(rmax, zval, xmin, dx, ibound)

    @printf("Nrmesh = %d\n", grid.Nrmesh)
    for i in 1:10
        @printf("%5d %18.10f\n", i, grid.r[i])
    end

    println("Pass here ...")
end

main()