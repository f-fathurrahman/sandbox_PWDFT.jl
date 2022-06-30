#---------------------------------------
abstract type AbstractXC end

struct Libxc <: AbstractXC
    x::Float64
end

struct XCDefault <: AbstractXC
    y::Float64
end
#---------------------------------------



# ----------------------------------------
abstract type AbstractPsPot end

struct PsPot_GTH <: AbstractPsPot
    x::Float64
end

struct PsPot_UPF <: AbstractPsPot
    y::Float64
end
# ----------------------------------------




struct SuperStruct{T1 <: AbstractXC, T2 <: AbstractPsPot}
    va::T1
    vb::Vector{T2}
end


my_str1 = SuperStruct(
    XCDefault(0.1),
    [PsPot_UPF(1.0), PsPot_UPF(2.0)]
)


my_str2 = SuperStruct(
    Libxc(9.9),
    [PsPot_GTH(9.9), PsPot_GTH(1.0)]
)

println(my_str1)

# The most general
function my_func1( m::SuperStruct )
    println("This should be: (Libxc,XCDefault) and (PsPot_GTH,PsPot_UPF)")
end

function my_func1( m::SuperStruct{T,PsPot_UPF} ) where T <: AbstractXC
    println("This should be: (Libxc,XCDefault) and (PsPot_UPF)")
end

# where T1 <: AbstractXC where T2 <: AbstractPsPot
# where {T1 <: AbstractXC, T2 <: AbstractPsPot}