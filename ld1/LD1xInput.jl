mutable struct LD1XInput
    Zval::Float64
    Zed::Float64 # redundant?
    Nspin::Int64
    Nwf::Int64
    nn::Vector{Int64}
    ll::Vector{Int64}
    oc::Vector{Float64}
    iswitch::Int64
    rel::Int64
    rcloc::Float64
    el::Vector{String}
    els::Vector{String}
    core_state::Vector{Bool}
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
    Nwfs =  6
    el = ["1S", "2S", "2P", "3S", "3P", "3D"] 
    els = ["3S", "3S", "3P", "3P", "3D", "3D"] 

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf
    @assert length(el) == Nwf
    @assert length(els) == Nwfs

    core_state = _calc_core_state(el, els)

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf,
        nn, ll, oc,
        iswitch, rel,
        rcloc,
        el, els, core_state
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
    els = ["5S", "5S", "5P", "5P", "4D", "4D"]

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf
    @assert length(el) == Nwf
    @assert length(els) == Nwfs

    core_state = _calc_core_state(el, els)

    # FIXME: define isw: spin index
    rel = 1
    iswitch = 1
    rcloc = 2.2

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf,
        nn, ll, oc,
        iswitch, rel,
        rcloc, 
        el, els, core_state
    )
end