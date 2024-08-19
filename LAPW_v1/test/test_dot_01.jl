using LinearAlgebra
N = 4
a = zeros(ComplexF64, N)
b = zeros(ComplexF64, N)

fill!(a, 1.0 + 2im)
fill!(b, 1.0 + 2im)
println(dot(a, b))

fill!(a, 1.0 - 2im)
fill!(b, 1.0 - 2im)
println(dot(a, b))
println(dot(b, a))
println(dot(conj(a), b))
println(dot(a, conj(b)))
