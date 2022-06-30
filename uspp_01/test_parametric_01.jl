abstract type AbstractT1 end
abstract type AbstractT2 end

struct T1A <: AbstractT1
    s::Symbol
    x::Float64
end

struct T1B <: AbstractT1
    s::Symbol
    y::Float64
end

struct T2A <: AbstractT2
    ss::Symbol
    xx::Float64
end

struct T2B <: AbstractT2
    ss::Symbol
    yy::Float64
end

struct SuperStruct{T1 <: AbstractT1, T2 <: AbstractT2}
    va::Vector{T1}
    vb::Vector{T2}
end


my_str1 = SuperStruct(
    [T1B(:b_name1, 1.1), T1B(:b_name2, 2.2)],
    [T2A(:a_name1, 2.3), T2A(:a_name2, 3.9)]
)


my_str2 = SuperStruct(
    [T1A(:b_name1, 1.1), T1A(:b_name2, 2.2)],
    [T2B(:a_name1, 2.3), T2B(:a_name2, 3.9)]
)

println(my_str1)

function my_func1( m::SuperStruct{T1A,T2B} )
    println("This should be T1A and T2B")
end


function my_func1( m::SuperStruct{T,T2B} ) where T <: AbstractT1
    println("This should be either (T1A or T1B) and T2B")
end

function my_func1( m::SuperStruct )
    println("This should be either (T1A or T1B) and (T1A or T2B)")
end