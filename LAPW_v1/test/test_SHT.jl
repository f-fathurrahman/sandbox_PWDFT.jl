using LAPWDFT

const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/LAPW/src/liblapwdft.so"

function get_elk_arrays(lmaxi::Int64, lmaxo::Int64)

    lmmaxi = (lmaxi + 1)^2
    lmmaxo = (lmaxo + 1)^2

    # Initialize, using custom subroutine
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

    ptr = cglobal( (:__m_sht_MOD_rbshto, LIBLAPW), Ptr{Float64} )
    tmp = zeros(Float64,lmmaxo*lmmaxo)
    ip = 1
    for j in 1:lmmaxo, i in 1:lmmaxo
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
        ip = ip + 1
    end
    rbshto = reshape(tmp, (lmmaxo,lmmaxo))

    return rbshti, rbshto
end


lmaxi = 1
lmaxo = 6
shtmat = SphericalHarmonicTransform(lmaxi, lmaxo)

rbshti, rbshto = get_elk_arrays(lmaxi, lmaxo)

# Should be zeros
display(shtmat.rbshti - rbshti); println()
display(shtmat.rbshto - rbshto); println()

