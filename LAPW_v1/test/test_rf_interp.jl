using Printf

include("../src/rf_interp.jl")
include("../src/spline.jl")

function my_func(x)
    return exp(-0.1*x)*sin(x)
end

function main()
    ni = 30
    xi = zeros(Float64,ni)
    fi = zeros(Float64,ni)

    no = 7
    xo = zeros(Float64,no)
    fo = zeros(Float64,no)

    L = 1.0
    dx = L/(ni-1)
    for i in 1:ni
        xi[i] = (i-1)*dx
        fi[i] = my_func(xi[i])
        @printf("%4d %18.10f %18.10f\n", i, xi[i], fi[i])
    end

    dx = L/(no-1)
    for i in 1:no
        xo[i] = (i-1)*dx
    end
    xo[end] = L-0.01

    rf_interp!(ni, xi, fi, no, xo, fo)

    println("Interpolation")
    for i in 1:no
        fxo = my_func(xo[i])
        Δ = abs(fxo - fo[i])
        @printf("%4d %18.10f %18.10f %18.10f %18.10e\n", i, xo[i], fo[i], fxo, Δ)
    end

end

@time main()
@time main()