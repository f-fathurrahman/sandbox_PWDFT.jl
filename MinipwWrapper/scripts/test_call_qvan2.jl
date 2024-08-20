using Printf

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )

ngm = unsafe_load(cglobal((:__gvect_MOD_ngm, LIBMINIPW), Int32)) |> Int64
println("ngm = ", ngm)

tpiba = unsafe_load(cglobal((:__cell_base_MOD_tpiba, LIBMINIPW), Float64))
println("tpiba = ", tpiba)

lmaxq = unsafe_load(cglobal((:__uspp_param_MOD_lmaxq, LIBMINIPW), Int32)) |> Int64
println("lmaxq = ", lmaxq)

# Read G-vectors
ptr = cglobal( (:__gvect_MOD_g, LIBMINIPW), Ptr{Float64} )
tmp = zeros(Float64, 3*ngm)
ip = 1
for j in 1:ngm, i in 1:3
    tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
    ip = ip + 1
end
G = reshape(tmp,3,ngm)

# Read G2
ptr = cglobal( (:__gvect_MOD_gg, LIBMINIPW), Ptr{Float64})
G2 = zeros(Float64, ngm)
for i in 1:ngm
    G2[i] = unsafe_load(unsafe_load(ptr,1),i)
end

# Calculate modulus of G-vectors, convert to atomic unit
Gmod = zeros(Float64, ngm)
for i in 1:ngm
    Gmod[i] = sqrt(G2[i])*tpiba
end

QfuncG = zeros(ComplexF64, ngm)
ylmk0 = zeros(Float64, ngm, lmaxq*lmaxq)


# CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
lmmax = lmaxq*lmaxq
ccall(
    (:ylmr2_, LIBMINIPW), Cvoid, 
    (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
    Int32(lmmax), Int32(ngm), G, G2, ylmk0
)


println(ylmk0[4,:])

# Input
isp = 2
ih = 1
jh = 1
# TODO: check sanity of the input

ccall(
    (:my_qvan2_, LIBMINIPW), Cvoid,
    (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32},
     Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}),
    Int32(ngm), Int32(ih), Int32(jh), Int32(isp),
    Gmod, QfuncG, ylmk0
)

for i in 1:5
    @printf("%5d [%18.10f,%18.10f]\n", i, real(QfuncG[i]), imag(QfuncG[i]))
end

println("sum QfuncG = ", sum(abs.(QfuncG)))