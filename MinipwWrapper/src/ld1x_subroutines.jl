# we assume that either pw.x or ld1.x is called
function call_ld1x_prepare_all(; filename=nothing)
    # XXX We also use global variable IS_PREPARED
    if MinipwWrapper.IS_PREPARED
        error("This function should be called only ONCE")
    end
    # If running from REPL filename must be given, otherwise
    # it will wait for stdin
    if isinteractive()
        if isnothing(filename)
            error("\n\n **** Please pass filename kwarg as input ****\n\n")
        end
    end
    #
    # If filename is given, redirect IOStream read from that file to stdin
    if !isnothing(filename)
        _call_func() = ccall( (:ld1x_prepare_all_, LIBMINIPW), Cvoid, () )
        file = open(filename, "r")
        redirect_stdin(_call_func, file)
        close(file)
    end
    MinipwWrapper.IS_PREPARED = true
    return
end

# actually calls driver_starting_potential
function call_ld1x_starting_potential()
    ccall( (:ld1x_driver_starting_potential_, LIBMINIPW), Cvoid, () )
end
