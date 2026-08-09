#XXX Some of these variables are not really input variables
mutable struct LD1XInput
    Zval::Float64
    Zed::Float64 # redundant?
    Nspin::Int64
    Nwf::Int64
    Nwfs::Int64
    nn::Vector{Int64}
    ll::Vector{Int64}
    oc::Vector{Float64}
    iswitch::Int64
    rel::Int64
    rcloc::Float64
    el::Vector{String}
    els::Vector{String}
    core_state::Vector{Bool}
    rcore::Float64
    nstoae::Vector{Int64}
    nns::Vector{Int64}
    lls::Vector{Int64}
    ocs::Vector{Float64}
    Enls::Vector{Float64}
    rcut::Vector{Float64}
    rcutus::Vector{Float64}
    isws::Vector{Float64}
    Nbeta::Int64
    fit_to_arbitrary_energy::Vector{Bool}
end

function _calc_core_state(el, els)
    Nwf = size(el, 1)
    core_state = ones(Bool, Nwf)
    for iwf in 1:Nwf
        if el[iwf] in els
            core_state[iwf] = false
        end
    end
    return core_state
end

function _calc_nstoae(el, els)
    Nwf = size(el, 1)
    Nwfs = size(els, 1)
    nstoae = zeros(Int64, Nwfs)
    for iwfs in 1:Nwfs
        nstoae[iwfs] = 0
        for iwf in 1:Nwf
            # XXX for fully relativistic also need to compare jj and jjs
            if els[iwfs] == el[iwf]
                nstoae[iwfs] = iwf
            end
        end
    end
    return nstoae
end


function _read_ps_config(str_config)
    lines = split(str_config, "\n", keepempty = false)
    Nwfs = parse(Int64, lines[1])
    els = Vector{String}(undef, Nwfs)
    nns = Vector{Int64}(undef, Nwfs)
    lls = Vector{Int64}(undef, Nwfs)
    ocs = Vector{Float64}(undef, Nwfs)
    Enls = Vector{Float64}(undef, Nwfs)
    rcut = Vector{Float64}(undef, Nwfs)
    rcutus = Vector{Float64}(undef, Nwfs)
    isws = Vector{Float64}(undef, Nwfs)
    iline = 2
    for iwfs in 1:Nwfs
        ll = split(lines[iline], " ", keepempty = false)
        els[iwfs] = ll[1]
        nns[iwfs] = parse(Int64, ll[2])
        lls[iwfs] = parse(Int64, ll[3])
        ocs[iwfs] = parse(Float64, ll[4])
        Enls[iwfs] = 0.5*parse(Float64, ll[5]) #XXX Convert to Ha
        rcut[iwfs] = parse(Float64, ll[6])
        rcutus[iwfs] = parse(Float64, ll[7])
        isws[iwfs] = parse(Float64, ll[8])
        # Next line
        iline += 1
    end
    return Nwfs, els, nns, lls, ocs, Enls, rcut, rcutus, isws
end




# Si uspp
function create_input_Si()

    Zval = 14.0
    Zed = Zval # no need to convert it to Ry, how about the sign?
    Nspin = 1
    Nwf = 6
    nn = [1, 2, 2, 3, 3, 3]
    ll = [0, 0, 1, 0, 1, 2]
    oc = [2.0, 2.0, 6.0, 2.0, 2.0, -1.0] # why -1 ?

    # FIXME: define isw: spin index
    rel = 1
    iswitch = 1
    rcloc = 1.9
    rcore = 1.3
    el = ["1S", "2S", "2P", "3S", "3P", "3D"] 

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf
    @assert length(el) == Nwf

    str_ps_config = """
    6
    3S  1  0  2.00  0.00  1.60  1.80  0.0
    3S  1  0  0.00  6.00  1.60  1.80  0.0
    3P  2  1  2.00  0.00  1.60  1.80  0.0
    3P  2  1  0.00  6.00  1.60  1.80  0.0
    3D  3  2  0.00  0.10  1.60  1.80  0.0
    3D  3  2  0.00  0.30  1.60  1.80  0.0
    """

    Nwfs, els, nns, lls, ocs, Enls, rcut, rcutus, isws = _read_ps_config(str_ps_config)
    
    @assert length(els) == Nwfs

    core_state = _calc_core_state(el, els)
    nstoae = _calc_nstoae(el, els)

    #nsloc = -1
    Nbeta = Nwfs # XXX This is only for lloc == -1

    SMALL_ENERGY = 1e-13
    fit_to_arbitrary_energy = zeros(Bool, Nwfs)
    for iwfs in 1:Nwfs
        if Enls[iwfs] >= SMALL_ENERGY # or equal to zero
            fit_to_arbitrary_energy[iwfs] = true
        end
    end

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf, Nwfs,
        nn, ll, oc,
        iswitch, rel,
        rcloc,
        el, els, core_state,
        rcore, nstoae,
        nns, lls, ocs, Enls, rcut, rcutus, isws,
        Nbeta,
        fit_to_arbitrary_energy
    )

end

function create_input_Pd()
    Zval = 46.0
    Zed = Zval
    Nspin = 1
    Nwf = 11
    nn = [1, 2, 2, 3, 3, 4, 4, 3, 5, 5, 4]
    ll = [0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2]
    oc = [2.0, 2.0, 6.0, 2.0, 6.0, 2.0, 6.0, 10.0, 1.0, 0.0, 9.0]
    Nwfs = 6
    el = ["1S", "2S", "2P", "3S", "3P", "4S", "4P", "3D", "5S", "5P", "4D"]

    # FIXME: define isw: spin index
    rel = 1
    iswitch = 1
    rcloc = 2.2
    rcore = 1.8

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf
    @assert length(el) == Nwf

    str_ps_config = """
    6
    5S 1 0 1.00  0.00  2.00  2.40  0.0
    5S 1 0 0.00  4.00  2.00  2.40  0.0
    5P 2 1 0.00  0.00  2.00  2.60  0.0
    5P 2 1 0.00  3.50  2.00  2.60  0.0
    4D 3 2 9.00  0.00  0.90  1.90  0.0
    4D 3 2 0.00  0.05  0.90  1.90  0.0
    """
    Nwfs, els, nns, lls, ocs, Enls, rcut, rcutus, isws = _read_ps_config(str_ps_config)

    @assert length(els) == Nwfs

    core_state = _calc_core_state(el, els)
    nstoae = _calc_nstoae(el, els)

    #nsloc = -1
    Nbeta = Nwfs # XXX This is only for lloc == -1

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf, Nwfs,
        nn, ll, oc,
        iswitch, rel,
        rcloc, 
        el, els, core_state,
        rcore, nstoae,
        nns, lls, ocs, Enls, rcut, rcutus, isws,
        Nbeta,
        fit_to_arbitrary_energy
    )
end

