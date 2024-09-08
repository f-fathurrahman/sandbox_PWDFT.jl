module ElkDFTWrapper_v632

using PWDFT: KPoints
using OffsetArrays
using Serialization: serialize

import ElkDFTWrapper
const PATH_ELKDFTWRAPPER = joinpath(pathof(ElkDFTWrapper), "..")

# This is currently hardcoded
const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/LAPW/src_v6.3.2/liblapwdft.so"

const ELK_DATA_DIR = "./TEMP_datadir"


include(joinpath(PATH_ELKDFTWRAPPER, "load_array.jl"))

include(joinpath(PATH_ELKDFTWRAPPER, "lattice.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "atomic_species.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "muffin_tins.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "gvectors.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "gkvectors.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "atoms.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "symmetry.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "density_pot_xc.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "sht.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "core_states.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "apwlo.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "info_apwlo.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "spin.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "kpoints.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "hamiltonian.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "energy.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "states.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "convergence.jl"))
include(joinpath(PATH_ELKDFTWRAPPER, "charge_moment_current.jl"))

include(joinpath(PATH_ELKDFTWRAPPER, "serialize_variables.jl"))
export serialize_variables

include("subroutines.jl")

#=
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

function export_pwdft_kpoints(; filename="pwdftjl_kpoints.dat")
    Nkpt = get_nkpt()
    ngridk = get_ngridk()
    bvec = get_bvec()
    vkc = get_vkc()[:,1:Nkpt]
    wk = get_wkpt()[1:Nkpt]
    kpoints = KPoints(
        Nkpt, (ngridk[1],ngridk[2],ngridk[3]),
        vkc, wk, bvec
    )
    @info "PWDFT.jl KPoints will be serialized to $filename"
    serialize(filename, kpoints)
end

# Put the workload that we want to investigate
function init_debug_calc()
    
    init_run()
    export_pwdft_kpoints()
    
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
=#

end
