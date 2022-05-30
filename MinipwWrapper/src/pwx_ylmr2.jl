function pwx_ylmr2!(
    lmax::Int64,
    r::Matrix{Float64},
    Ylm::Matrix{Float64}
)

    @assert size(r,1) == 3
    N = size(r,2)
    lmmax = (lmax+1)^2
    Ylm = zeros(Float64, N, lmmax)
    rr = zeros(Float64, N)
    for i in 1:N
        rr[i] = r[1,i]^2 + r[2,i]^2 + r[3,i]^2
    end

    ccall(
        (:ylmr2_, LIBMINIPW), Cvoid, 
        (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        Int32(lmmax), Int32(N), r, rr, Ylm
    )

    return
end
