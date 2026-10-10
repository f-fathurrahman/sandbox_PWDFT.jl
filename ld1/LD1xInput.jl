
# zed: 


#XXX Some of these variables are not really input variables
mutable struct LD1XInput
    Zval::Float64  # the ionic pseudo charge
    Zed::Float64 # the ionic charge, XXX: redundant? the same as Zed?
    Nspin::Int64
    Nwf::Int64
    Nwf_ps::Int64
    nn::Vector{Int64}
    ll::Vector{Int64}
    oc::Vector{Float64}
    iswitch::Int64
    rel::Int64
    rcloc::Float64
    state_label::Vector{String}
    state_label_ps::Vector{String}
    core_state::Vector{Bool}
    rcore::Float64
    idx_ps_to_ae::Vector{Int64}
    nn_ps::Vector{Int64}
    ll_ps::Vector{Int64}
    oc_ps::Vector{Float64}
    Enl_ps::Vector{Float64}
    rcut::Vector{Float64}
    rcutus::Vector{Float64}
    isws::Vector{Float64}
    Nbeta::Int64
    fit_to_arbitrary_energy::Vector{Bool}
end

function _calc_core_state(state_label, state_label_ps)
    Nwf = size(state_label, 1)
    core_state = ones(Bool, Nwf)
    for iwf in 1:Nwf
        if state_label[iwf] in state_label_ps
            core_state[iwf] = false
        end
    end
    return core_state
end

function _calc_nstoae(state_label, state_label_ps)
    Nwf = size(state_label, 1)
    Nwf_ps = size(state_label_ps, 1)
    idx_ps_to_ae = zeros(Int64, Nwf_ps)
    for iwfs in 1:Nwf_ps
        idx_ps_to_ae[iwfs] = 0
        for iwf in 1:Nwf
            # XXX for fully relativistic also need to compare jj and jjs
            if state_label_ps[iwfs] == state_label[iwf]
                idx_ps_to_ae[iwfs] = iwf
            end
        end
    end
    return idx_ps_to_ae
end


function _read_ps_config(str_config)
    lines = split(str_config, "\n", keepempty = false)
    Nwf_ps = parse(Int64, lines[1])
    state_label_ps = Vector{String}(undef, Nwf_ps)
    nn_ps = Vector{Int64}(undef, Nwf_ps)
    ll_ps = Vector{Int64}(undef, Nwf_ps)
    oc_ps = Vector{Float64}(undef, Nwf_ps)
    Enl_ps = Vector{Float64}(undef, Nwf_ps)
    rcut = Vector{Float64}(undef, Nwf_ps)
    rcutus = Vector{Float64}(undef, Nwf_ps)
    isws = Vector{Float64}(undef, Nwf_ps)
    iline = 2
    for iwfs in 1:Nwf_ps
        ll = split(lines[iline], " ", keepempty = false)
        state_label_ps[iwfs] = ll[1]
        nn_ps[iwfs] = parse(Int64, ll[2])
        ll_ps[iwfs] = parse(Int64, ll[3])
        oc_ps[iwfs] = parse(Float64, ll[4])
        Enl_ps[iwfs] = 0.5*parse(Float64, ll[5]) #XXX Convert to Ha
        rcut[iwfs] = parse(Float64, ll[6])
        rcutus[iwfs] = parse(Float64, ll[7])
        isws[iwfs] = parse(Float64, ll[8])
        # Next line
        iline += 1
    end
    return Nwf_ps, state_label_ps, nn_ps, ll_ps, oc_ps, Enl_ps, rcut, rcutus, isws
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
    state_label = ["1S", "2S", "2P", "3S", "3P", "3D"] 

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf
    @assert length(state_label) == Nwf

    str_ps_config = """
    6
    3S  1  0  2.00  0.00  1.60  1.80  0.0
    3S  1  0  0.00  6.00  1.60  1.80  0.0
    3P  2  1  2.00  0.00  1.60  1.80  0.0
    3P  2  1  0.00  6.00  1.60  1.80  0.0
    3D  3  2  0.00  0.10  1.60  1.80  0.0
    3D  3  2  0.00  0.30  1.60  1.80  0.0
    """

    Nwf_ps, state_label_ps, nn_ps, ll_ps, oc_ps, Enl_ps, rcut, rcutus, isws = _read_ps_config(str_ps_config)
    
    @assert length(state_label_ps) == Nwf_ps

    core_state = _calc_core_state(state_label, state_label_ps)
    idx_ps_to_ae = _calc_nstoae(state_label, state_label_ps)

    #nsloc = -1
    Nbeta = Nwf_ps # XXX This is only for lloc == -1

    SMALL_ENERGY = 1e-13
    fit_to_arbitrary_energy = zeros(Bool, Nwf_ps)
    for iwfs in 1:Nwf_ps
        if Enl_ps[iwfs] >= SMALL_ENERGY # or equal to zero
            fit_to_arbitrary_energy[iwfs] = true
        end
    end

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf, Nwf_ps,
        nn, ll, oc,
        iswitch, rel,
        rcloc,
        state_label, state_label_ps, core_state,
        rcore, idx_ps_to_ae,
        nn_ps, ll_ps, oc_ps, Enl_ps, rcut, rcutus, isws,
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
    state_label = ["1S", "2S", "2P", "3S", "3P", "4S", "4P", "3D", "5S", "5P", "4D"]

    # FIXME: define isw: spin index
    rel = 1
    iswitch = 1
    rcloc = 2.2
    rcore = 1.8

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf
    @assert length(state_label) == Nwf

    str_ps_config = """
    6
    5S 1 0 1.00  0.00  2.00  2.40  0.0
    5S 1 0 0.00  4.00  2.00  2.40  0.0
    5P 2 1 0.00  0.00  2.00  2.60  0.0
    5P 2 1 0.00  3.50  2.00  2.60  0.0
    4D 3 2 9.00  0.00  0.90  1.90  0.0
    4D 3 2 0.00  0.05  0.90  1.90  0.0
    """
    Nwf_ps, state_label_ps, nn_ps, ll_ps, oc_ps, Enl_ps, rcut, rcutus, isws = _read_ps_config(str_ps_config)

    @assert length(state_label_ps) == Nwf_ps

    core_state = _calc_core_state(state_label, state_label_ps)
    idx_ps_to_ae = _calc_nstoae(state_label, state_label_ps)

    #nsloc = -1
    Nbeta = Nwf_ps # XXX This is only for lloc == -1

    SMALL_ENERGY = 1e-13
    fit_to_arbitrary_energy = zeros(Bool, Nwf_ps)
    for iwfs in 1:Nwf_ps
        if Enl_ps[iwfs] >= SMALL_ENERGY # or equal to zero
            fit_to_arbitrary_energy[iwfs] = true
        end
    end

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf, Nwf_ps,
        nn, ll, oc,
        iswitch, rel,
        rcloc, 
        state_label, state_label_ps, core_state,
        rcore, idx_ps_to_ae,
        nn_ps, ll_ps, oc_ps, Enl_ps, rcut, rcutus, isws,
        Nbeta,
        fit_to_arbitrary_energy
    )
end

