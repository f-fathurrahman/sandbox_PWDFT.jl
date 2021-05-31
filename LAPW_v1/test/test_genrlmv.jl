using Printf
using LAPWDFT

function main()
    v = [1.0, 1.1, 1.2]
    lmax = 2
    lmmax = (lmax+1)^2
    rlm = zeros(Float64,lmmax)
    genrlmv!(lmax, v, rlm)
    for i in 1:lmmax
        @printf(" %4d%18.10f\n", i, rlm[i])
    end
end

main()