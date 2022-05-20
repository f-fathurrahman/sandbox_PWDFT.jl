module MinipwWrapper

using OffsetArrays

# This is currently hardcoded
const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

function pwx_prepare_all()
    ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )
    return
end



end # module
