function call_initialize_elk()
    ccall( (:initialize_elk_, LIBLAPW), Cvoid, () )
    return
end

function call_init0()
    ccall( (:init0_, LIBLAPW), Cvoid, () )
    return
end

function call_init1()
    ccall( (:init1_, LIBLAPW), Cvoid, () )
    return
end

function call_rhoinit()
    ccall( (:rhoinit_, LIBLAPW), Cvoid, () )
    return
end

function call_maginit()
    ccall( (:maginit_, LIBLAPW), Cvoid, () )
    return
end

function call_potks(; txc=true)
    ccall( (:potks_, LIBLAPW), Cvoid, (Ref{Bool},), txc )
    return
end

# This is the original gndstate
function call_gndstate()
    ccall( (:gndstate_, LIBLAPW), Cvoid, () )
    return
end