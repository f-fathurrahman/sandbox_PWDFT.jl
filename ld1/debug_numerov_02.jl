# Based on: http://staff.ustc.edu.cn/~zqj/posts/Numerov-Algorithm/

using Printf

function radial_wfc_numerov_at0(E, r0; n=1, l=0, Z=1, du=0.001)
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

    # Backward integration, start from practical Inf
    ur[end] = 0.0
    ur[end-1] = du # FIXME: set to a small quantity?

    dr  = r0[2] - r0[1]
    h12 = dr^2/12.0

    gn = zeros(Float64, Nr)
    fn = zeros(Float64, Nr)
    for i in 1:Nr
        gn[i] = E + 2*Z / r0[i] - l*(l+1)/r0[i]^2
        fn[i] = 1.0 + h12 * gn[i]
    end

    for ii in range(Nr-2, stop=1, step=-1)
        ur[ii] = (12 - 10*fn[ii+1]) * ur[ii+1] - ur[ii+2] * fn[ii+2]
        ur[ii] = ur[ii]/fn[ii]
    end

    # normalization, using simple rectangle rule
    ss = sum(ur.^2) * dr
    @views ur[:] .= ur[:]/sqrt(ss)

    # now extrapolate the wavefunction to u(0)
    # recall that radial grid does not include 0.
    # the first derivative equals at ur[0]
    #
    # (u0 - ur[1]) / (0 - r0[1]) = (ur[2] - ur[1]) / (r0[2] - r0[1])
    # r0[2] - r0[1] = dr
    #
    u0 = ur[1] + (ur[2] - ur[1])/dr * (0 - r0[1])

    return u0
end


function my_bisection(f, x1, x2; TOL=1.0e-9, verbose=false, NiterMax=nothing)

    if verbose
        println("")
        println("Searching root with bisection method:")
        @printf("Interval = (%18.10f,%18.10f)\n", x1,x2)
        @printf("TOL = %e\n", TOL)
        println("")
    end

    f1 = f(x1)
    if abs(f1) <= TOL
        return x1, 0.0
    end

    f2 = f(x2)
    if abs(f2) <= TOL
        return x2, 0.0
    end

    if f1*f2 > 0.0
        error("Root is not bracketed")
    end

    # No NiterMax is provided
    # We calculate the default value here.
    if NiterMax == nothing
        NiterMax = round(Int64, ceil(log10(abs(x2-x1)/TOL))/log10(2.0)) + 10
        # extra 10 iterations
    end
    println("NiterMax = ", NiterMax)

    # For the purpose of calculating relative error
    x3 = 0.0
    x3_old = 0.0

    for i in range(1,NiterMax+1)

        x3_old = x3
        x3 = 0.5*(x1 + x2)
        f3 = f(x3)

        if verbose
            @printf("bisection: %5d %18.10f %15.5e\n", i, x3, abs(f3))
        end

        if abs(f3) <= TOL
            if verbose
                println("")
                @printf("bisection is converged in %d iterations\n", i)
            end
            # return the result
            return x3, abs(x3 - x3_old)
        end


        if f2*f3 < 0.0
            # sign of f2 and f3 is different
            # root is in [x2,x3]
            # change the interval bound of x1 to x3
            x1 = x3
            f1 = f3
        else
            # sign of f1 and f3 is different
            # root is in [x1,x3]
            # change the interval bound of x2 to x3
            x2 = x3
            f2 = f3
        end
    end


    # No root is found after NiterMax iterations
    if verbose
        println("No root is found")
    end
    return nothing, nothing
end


function main()
    r0 = range(1e-7, stop=30.0, length=3000)

    Z = 2
    l = 0
    n = 2

    E_lower = -0.3
    E_upper = E_lower
    dE = 0.10
    u1 = radial_wfc_numerov_at0(E_lower, r0, n=n, l=l, Z=Z)

    # Search for E_upper
    while true
        E_upper += dE # increase
        u2 = radial_wfc_numerov_at0(E_upper, r0, n=n, l=l, Z=Z)
        if u1 * u2 < 0
            break
        end
    end

    my_func(E) = radial_wfc_numerov_at0(E, r0, n=n, l=l, Z=Z, du=0.001)
    Eroot, _ = my_bisection(my_func, E_lower, E_upper, verbose=true)
    println(Eroot)
end


main()

