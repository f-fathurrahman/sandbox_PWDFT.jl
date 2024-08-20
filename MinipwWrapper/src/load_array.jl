# Load allocatable array of eltype T, 
# dimensions given in `shape` and gfortran `symbol`
function _load_allocatable_array(symbol::Symbol, T::DataType, shape)
    # Special case for Int64
    if T == Int64
        TT = Int32
    else
        TT = T
    end
    ptr = cglobal( (symbol, LIBMINIPW), Ptr{TT} )
    tmp = zeros(T, prod(shape))
    for ip in 1:prod(shape)
        tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
    end
    return reshape(tmp, shape)
end


# Not tested for multidimensional automatic array
function _load_automatic_array(symbol::Symbol, T::DataType, shape)
    # Special case for Int64
    if T == Int64
        TT = Int32
    else
        TT = T
    end
    ptr = cglobal( (symbol, LIBMINIPW), TT )
    tmp = zeros(T, prod(shape))
    for ip in 1:prod(shape)
        tmp[ip] = unsafe_load(ptr,ip)
    end
    return reshape(tmp, shape)
end