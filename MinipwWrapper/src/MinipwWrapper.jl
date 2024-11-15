module MinipwWrapper

using OffsetArrays

# This is currently hardcoded
const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

IS_PREPARED = false

# !!!! This function should be called only once
function call_prepare_all(; filename=nothing)
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
        _call_func() = ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )
        file = open(filename, "r")
        redirect_stdin(_call_func, file)
        close(file)
    end
    MinipwWrapper.IS_PREPARED = true
    return
end

include("load_array.jl")
include("subroutines.jl")
include("variables.jl")

include("ld1x_subroutines.jl")
include("ld1x_variables.jl")

# Do a single-point calculation, also calculate forces and stress
function run_single_point(; filename=nothing)
    #
    call_prepare_all(filename=filename)
    call_info_upf()
    call_my_electrons()
    call_my_forces()
    #
    stress = zeros(Float64,3,3)
    ccall( (:my_stress_, LIBMINIPW), Cvoid, (Ptr{Float64},), stress )
    #
    # Print out the stress tensor in original units for easy comparison
    println("\nStress obtained from my_stress\n")
    display(stress); println()
    #
    return
end



end # module
