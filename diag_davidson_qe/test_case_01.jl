using LinearAlgebra

const LIBDAVIDSON = "./libdavidson.so"

N = 20

H = rand(ComplexF64,N,N)
H = 0.5*(H + H')
H = H + N*I # Make diagonally dominant

S = rand(ComplexF64,N,N)
S = 0.5*(S + S')
S = S + N*I # Make diagonally dominant

res = eigen(H, S)
#println(res.values')

npw = N
npwx = N
nvec = 3
nvecx = 2*nvec
evc = randn(ComplexF64, npwx, nvec)
U = inv(sqrt(evc'*S*evc))
evc = evc*U
ethr = 1e-8
uspp = true
e = zeros(Float64, nvec)
btype = ones(Int32, nvec)
notcnv = Int32(0)
lrot = false
dav_iter = Int32(0)

ccall( (:my_cegterg_v61_, LIBDAVIDSON), Cvoid,
    (
        Ptr{ComplexF64}, Ptr{ComplexF64},
        Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ref{Float64},
        Ref{Bool}, Ptr{Float64}, Ptr{Int32}, Ref{Int32}, Ref{Bool}, Ref{Int32}
    ),
    H, S, 
    Int32(npw), Int32(npwx), Int32(nvec), Int32(nvecx), evc, ethr,
    uspp, e, btype, notcnv, lrot, dav_iter
)

println("res.values = ", res.values[1:nvec]')
println("e = ", e')

#SUBROUTINE my_cegterg_v61( &
#  H, S, & 
#  npw, npwx, nvec, nvecx, evc, ethr, &
#  uspp, e, btype, notcnv, lrot, dav_iter )
