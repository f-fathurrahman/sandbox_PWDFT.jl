using Printf

include("../src/sphcover.jl")
include("../src/sctovec.jl")

function main()
    n = 4
    tp = zeros(2,n)
    sphcover!(n, tp)
    println("θ and ϕ after sphcover for n = ", n)
    for i in 1:n
        @printf(" %8d%18.10f%18.10f\n", i, tp[1,i], tp[2,i])
    end

    v = zeros(3,n)
    sctovec!(n, tp, v)
    println("After sctovec for = ", n)
    for i in 1:n
        @printf(" %8d%18.10f%18.10f%18.10f\n", i, v[1,i], v[2,i], v[3,i])
    end
end

main()