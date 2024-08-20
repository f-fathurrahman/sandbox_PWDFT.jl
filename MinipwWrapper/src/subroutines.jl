function call_info_upf()
    ccall( (:info_upf_, LIBMINIPW), Cvoid, () )
    return
end

function call_my_electrons()
    ccall( (:my_electrons_, LIBMINIPW), Cvoid, () )
    return
end

function call_my_forces()
    ccall( (:my_forces_, LIBMINIPW), Cvoid, () )
    return
end

include("pwx_ylmr2.jl")
