using Printf

include("sphcover.jl")

function main()
    n = 10
    tp = zeros(2,10)
    sphcover!(n, tp)
    for i in 1:n
        @printf("%8d%18.10f%18.10f\n", i, tp[1,i], tp[2,i])
    end
end

main()