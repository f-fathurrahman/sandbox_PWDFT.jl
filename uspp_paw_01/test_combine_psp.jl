abstract type AbstractPspot end

struct PsPotFormatA <: AbstractPspot
end

struct PsPotFormatB <: AbstractPspot
end

struct ContainerPsPots{T<:AbstractPspot}
    pspots::Vector{T}
end

function main()
    pspA = PsPotFormatA()
    pspB = PsPotFormatB()
    pspots = ContainerPsPots([pspA, pspB])
    println("Pass here ...")
end

main()