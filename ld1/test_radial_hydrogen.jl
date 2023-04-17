using Printf
import PWDFT

include("RadialGrid.jl")

R10(r; Z=1, a0=1) = 2*(Z/a0)*exp(-Z*r/a0)

#function test_integ()

    zval = 1.0 # Depend on the atom
    rmax = 100.0
    xmin = -7.0 # iswitch = 1
    dx = 0.008
    ibound = false
    grid = RadialGrid(rmax, zval, xmin, dx, ibound)

    fr = zeros(Float64, grid.Nrmesh)
    for i in 1:grid.Nrmesh
        fr[i] = R10(grid.r[i])
    end

    res = PWDFT.integ_simpson(grid.Nrmesh, fr .* fr .* grid.r2, grid.rab)
    println("res = ", res)
#end
#test_integ()