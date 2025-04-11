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
end

function create_input_Si()
    Zval = 14.0
    Zed = Zval # no need to convert it to Ry, how about the sign?
    Nspin = 1
    Nwf = 5
    nn = [1, 2, 2, 3, 3]
    ll = [0, 0, 1, 0, 1] 
    oc = [2.0, 2.0, 6.0, 2.0, 2.0]
    # FIXME: define isw: spin index
    rel = 1
    iswitch = 1

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf,
        nn, ll, oc,
        iswitch, rel
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

    @assert length(nn) == Nwf
    @assert length(ll) == Nwf
    @assert length(oc) == Nwf

    # FIXME: define isw: spin index
    rel = 1
    iswitch = 1

    return LD1XInput(
        Zval, Zed,
        Nspin, Nwf,
        nn, ll, oc,
        iswitch, rel
    )
end