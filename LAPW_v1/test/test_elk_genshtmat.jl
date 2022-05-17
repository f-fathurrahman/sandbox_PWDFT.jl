const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/src_LAPW/liblapwdft.so"

lmaxi = 1
lmaxo = 6
ccall(
    (:my_genshtmat_, LIBLAPW), Cvoid, (Ref{Int32}, Ref{Int32}),
    Int32(lmaxi), Int32(lmaxo)
)

ptr = cglobal( (:__m_sht_MOD_rbshti, LIBLAPW), Ptr{Float64} )
lmmaxi = (lmaxi + 1)^2
tmp = zeros(Float64,lmmaxi*lmmaxi)
ip = 1
for j in 1:lmmaxi, i in 1:lmmaxi
    tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
    ip = ip + 1
end
rbshti = reshape(tmp, (lmmaxi,lmmaxi))
display(rbshti); println()


ptr = cglobal( (:__m_sht_MOD_rbshto, LIBLAPW), Ptr{Float64} )
lmmaxo = (lmaxo + 1)^2
tmp = zeros(Float64,lmmaxo*lmmaxo)
ip = 1
for j in 1:lmmaxo, i in 1:lmmaxo
    tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
    ip = ip + 1
end
rbshto = reshape(tmp, (lmmaxo,lmmaxo))
display(rbshto); println()


