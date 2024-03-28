# rlm: Rlm (real spherical harmonics) of size (lmax+1)^2
function genrlmv!(lmax, v, rlm)

    @assert lmax >= 0
    @assert lmax <= 50

    sqtwo = 1.4142135623730950488

    # generate complex spherical harmonics
    ylm = zeros(ComplexF64, (lmax+1)^2)
    genylmv!(lmax, v, ylm)
    
    # convert to real spherical harmonics
    lm = 0
    for l in 0:lmax
        for m in -l:-1
            lm = lm + 1
            rlm[lm] = sqtwo*imag(ylm[lm])
        end
        lm = lm + 1
        rlm[lm] = real(ylm[lm])
        for m in 1:l
            lm = lm + 1
            rlm[lm] = sqtwo*real(ylm[lm])
        end
    end
    return
end