# Based on: http://staff.ustc.edu.cn/~zqj/posts/Numerov-Algorithm/

function radial_wfc_numerov(r0; n=1, l=0, Z=1, du=0.001)
    # Numerov algorithm
    #
    #              [12 - 10f(n)]*y(n) - y(n-1)*f(n-1)
    #    y(n+1) = ------------------------------------
    #                           f(n+1)
    #
    # where
    #    
    #    f(n) = 1 + (h**2 / 12)*g(n)
    #    g(n) = [E + (2*Z / x) - l*(l+1) / x**2]

    # r0 is radial grid (linear)
    Nr = size(r0, 1)
    ur = zeros(Float64, Nr)
    fn = zeros(Float64, Nr)

    E = -Z^2/n^2
    println("E = ", E)

    # Backward integration, start from practical Inf
    ur[end] = 0.0
    ur[end-1] = du # FIXME: set to a small quantity?

    dr  = r0[2] - r0[1]
    h12 = dr^2/12.0

    gn = zeros(Float64, Nr)
    fn = zeros(Float64, Nr)
    for i in 1:Nr
        Veff = -2*Z/r0[i] + l*(l+1)/r0[i]^2
        gn[i] = E - Veff
        fn[i] = 1.0 + h12 * gn[i]
    end

    for ii in range(Nr-2, stop=1, step=-1)
        ur[ii] = (12 - 10*fn[ii+1]) * ur[ii+1] - ur[ii+2] * fn[ii+2]
        ur[ii] = ur[ii]/fn[ii]
    end

    # normalization
    ss = sum(ur.^2) * dr
    println("ss = ", ss)
    @views ur[:] .= ur[:]/sqrt(ss)

    return ur
end

import PyPlot
const plt = PyPlot

function main()
    r0 = range(1e-10, stop=20, length=500) |> collect

    Z = 2
    l = 0
    n = 2

    # backward integration from u(\infty) with Numerov method
    ur3 = radial_wfc_numerov(r0, n=n, l=l)
    
    plt.clf()
    plt.plot(r0, ur3)
    plt.grid(true)

    return
end

main()
