using Printf
using OffsetArrays

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

include("Ylm_real_qe.jl")
include("_generate_lm_indices_qe.jl")

Ng = 1
g = zeros(3,Ng)
g[:,1] = [0.1, 1.0, 2.0]
gg = [sum( g[:,1].^2) ]

lmax = 5
lmmax = (lmax+1)^2

idxlm, idxil, idxim = _generate_lm_indices_qe(lmax)
    
Ylm = zeros(1,lmmax)
@views Ylm_real_qe!(lmax, g[:], Ylm[1,:])

Ylm_ref = zeros(1,lmmax)
ccall(
    (:ylmr2_, LIBMINIPW), Cvoid, 
    (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
    Int32(lmmax), Int32(Ng), g, gg, Ylm_ref
)


# Using the "natural" loop (sum over l (sum over m))
#for l in 0:lmax
#    println()
#    for m in -l:l
#        lm = idxlm[l,m]
#        Δ = Ylm[1,lm] - Ylm_ref[1,lm]
#        @printf("l=%3d m=%3d lm=%3d Ylm = %18.10f %18.10f %18.10f\n",
#            l, m, lm, Ylm[1,lm], Ylm_ref[1,lm], Δ)
#    end
#end

println()
for lm in 1:lmmax
    l = idxil[lm]
    m = idxim[lm]
    Δ = Ylm[1,lm] - Ylm_ref[1,lm]
    @printf("l=%3d m=%3d lm=%3d Ylm = %18.10f %18.10f %18.10f\n",
            l, m, lm, Ylm[1,lm], Ylm_ref[1,lm], Δ)
end