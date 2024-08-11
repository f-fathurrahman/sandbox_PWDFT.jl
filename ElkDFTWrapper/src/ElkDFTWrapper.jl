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
include("gkvectors.jl")
include("atoms.jl")
include("symmetry.jl")
include("density_pot_xc.jl")
include("sht.jl")
include("core_states.jl")
include("apwlo.jl")
include("info_apwlo.jl")
include("spin.jl")
include("kpoints.jl")
include("hamiltonian.jl")
include("energy.jl")
include("states.jl")
include("convergence.jl")

include("serialize_variables.jl")
export serialize_variables


include("subroutines.jl")


# Need to call init_run separately
function gndstate()
    call_rhoinit()
    call_maginit()
    call_potks(txc=true)
    call_genvsig()
    call_my_gndstate_setup_mixing()

    # SCF loop
    # XXX some steps are removed
    Etot = get_engytot()
    Etot_old = Etot
    dE = NaN
    for iterSCF in 1:20
        call_my_gndstate_increment_iscl()
        call_gencore()
        call_linengy()
        call_genapwlofr()
        call_genevfsv()
        call_occupy()
        call_rhomag()
        call_potks(txc=true)
        call_my_gndstate_do_mixing()
        call_genvsig()
        call_energy()
        if iterSCF > 0
            Etot = get_engytot()
            dE = abs(Etot - Etot_old)
        end
        @info "iscl = $(get_iscl()) finished"
        @info "dE = $dE"
        Etot_old = Etot
    end
end


# Put the workload that we want to investigate
function init_debug_calc()
    
    init_run()
    
    # Only call after init_run
    #call_my_gndstate(1)
    # Otherwise call other debug subroutines below

    call_rhoinit()
    
    # potks or potks_no_symm
    call_potks()
    #call_potks_no_symm() # only for debugging
    
    call_genvsig()

    call_gencore()

    call_linengy()

    call_genapwfr()
    call_genlofr()
    call_olprad()
    call_hmlrad()

    call_genevfsv()

end
export init_debug_calc

end
