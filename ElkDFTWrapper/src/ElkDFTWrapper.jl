module ElkDFTWrapper

using OffsetArrays
using Serialization: serialize

# This is currently hardcoded
const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/LAPW/src/liblapwdft.so"

const ELK_DATA_DIR = "./TEMP_datadir"

function init_run()
    ccall( (:read_input_, LIBLAPW), Cvoid, () )
    ccall( (:init0_, LIBLAPW), Cvoid, () )
    ccall( (:init1_, LIBLAPW), Cvoid, () )
    ccall( (:info_gvectors_, LIBLAPW), Cvoid, () )
    ccall( (:info_muffin_tins_, LIBLAPW), Cvoid, () )
    ccall( (:writesym_, LIBLAPW), Cvoid, () )

    # Create the directory to save variables here
    isdir(ELK_DATA_DIR) || mkdir(ELK_DATA_DIR)
    return
end

function elk_rhoinit()
    ccall( (:rhoinit_, LIBLAPW), Cvoid, () )
    return
end

include("atomic_species.jl")
include("muffin_tins.jl")

include("atom.jl")
export elk_solve_atom!

include("density_pot_xc.jl")

include("write_elk_variables.jl")
export write_elk_variables

end
