
# From AERadialGridDescriptor in GPAW
struct AERadialGrid
    npts::Int64
    npts_spline::Int64
    a::Float64
    b::Float64
    r::Vector{Float64}
    dr::Vector{Float64}
    d2gdr2::Vector{Float64} # rename?
end

function AERadialGrid(
    a::Float64,
    b::Float64;
    npts=1000,
    npts_spline=25
)
    r = zeros(Float64, npts)
    dr = zeros(Float64, npts)
    d2gdr2 = zeros(Float64, npts)
    for i in 1:npts
        g = i - 1 # offset by 1
        r[i] = a*g/(1 - b*g)
        dr[i] = (b*r[i] + a)^2 / a
        d2gdr2[i] = -2 * a * b / (b * r[i] + a)^3
    end

    return AERadialGrid(
        npts, npts_spline,
        a, b, r, dr, d2gdr2
    )
end

import Base: print
function print(io::IO, obj::AERadialGrid)
    print(io, "AERadialGrid instance\n")
    print(io, "a=$(obj.a) b=$(obj.b) npts=$(obj.npts)")
end
