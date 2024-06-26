module ElkDFTWrapper

using OffsetArrays
using Serialization: serialize

# This is currently hardcoded
const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/LAPW/src/liblapwdft.so"

const ELK_DATA_DIR = "./TEMP_datadir"


include("load_array.jl")

"""
This should be run first. This will call some Elk subroutines
such as read_input, init0, init1 (and possibly others) such
that the global variables are allocated and (some of them) initilized.
This will also create ELK_DATA_DIR that will be used to save the
global variables.

NOTE: probably this must be called only once to avoid double allocations.
"""
function init_run()
    ccall( (:read_input_, LIBLAPW), Cvoid, () )
    ccall( (:init0_, LIBLAPW), Cvoid, () )
    ccall( (:init1_, LIBLAPW), Cvoid, () )

    # Create the directory to save variables here
    isdir(ELK_DATA_DIR) || mkdir(ELK_DATA_DIR)
    return
end

include("lattice.jl")
include("atomic_species.jl")
include("muffin_tins.jl")
include("gvectors.jl")
include("atoms.jl")
include("symmetry.jl")
include("density_pot_xc.jl")
include("sht.jl")

include("serialize_variables.jl")
export serialize_variables


include("subroutines.jl")


# Put the workload that we want to investigate
function init_debug_calc()
    init_run()
    call_rhoinit()
    call_potks()
end
export init_debug_calc

end
