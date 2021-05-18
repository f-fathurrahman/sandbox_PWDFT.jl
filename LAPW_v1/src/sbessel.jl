function sbessel!(
    lmax::Int64, x::Float64, jl
)
    for l in 0:lmax
        jl[l] = sphericalbesselj(l, x)
    end
    return
end