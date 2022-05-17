const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/src_LAPW/liblapwdft.so"

ccall( (:read_input_, LIBLAPW), Cvoid, () )
ccall( (:init0_, LIBLAPW), Cvoid, () )
ccall( (:init1_, LIBLAPW), Cvoid, () )

ptr = cglobal( (:__m_sht_MOD_rbshti, LIBLAPW), Ptr{Float64} )
lmmaxi = 4
tmp = zeros(Float64,lmmaxi*lmmaxi)
ip = 1
for j in 1:lmmaxi, i in 1:lmmaxi
    tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
    ip = ip + 1
end
rbshti = reshape(tmp, (lmmaxi,lmmaxi))
display(rbshti); println()