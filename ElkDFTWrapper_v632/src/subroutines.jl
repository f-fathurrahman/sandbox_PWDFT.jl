function call_initialize_elk()
    ccall( (:initialize_elk_, LIBLAPW), Cvoid, () )
    return
end

# This is the original gndstate
function call_gndstate()
    ccall( (:gndstate_, LIBLAPW), Cvoid, () )
    return
end