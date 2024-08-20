module MinipwWrapper

using OffsetArrays

# This is currently hardcoded
const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

function pwx_prepare_all(; filename=nothing)
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
        redirect_stdio(stdin=open(filename, "r"))
    end
    ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )
    return
end

function other_func()
    @info "Testing again ..."
end


end # module
