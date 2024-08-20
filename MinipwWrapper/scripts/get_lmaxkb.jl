using Serialization: serialize

const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )


function pwx_get_lmax()

    # Those are parameters (HARDCODED)
    lmaxx  = 3         # max non local angular momentum (l=0 to lmaxx)      
    lqmax = 2*lmaxx+1  # max number of angular momenta of Q

    lmaxkb = unsafe_load(cglobal((:__uspp_param_MOD_lmaxkb, LIBMINIPW), Int32)) |> Int64

    println("lmmaxx = ", lmaxx)
    println("lqmax = ", lqmax)
    println("lmaxkb = ", lmaxkb)

    # maximum number of combined angular momentum
    nlx = (lmaxx + 1)^2
    # maximum magnetic angular momentum of Q
    mx = 2*lqmax - 1

    println("nlx = ", nlx)
    println("mx  = ", mx)
end

pwx_get_lmax()
