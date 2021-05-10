using Printf
using SpecialFunctions: besselj

sbessel(l, x) = sqrt(pi/(2*x))*besselj(l+1/2, x)

function main()
    lmax = 4
    x = 1.2
    println("x = ", x)
    for l in 0:lmax
        res = sbessel(l, x)
        @printf(" %4d%18.10f\n", l, res)
    end
end

main()