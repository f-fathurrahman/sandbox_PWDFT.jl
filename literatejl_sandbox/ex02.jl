#=
# Introduction

In this introduction, we will see some Julia codes
=#

#=
First, we will import some packages
=#
using LinearAlgebra

#=
This is some simple calculations
=#
A = randn(3,3)
A = 0.5*(A' + A) # symmetrize
B = randn(3,3)
#Calculate determinant of A
det(A)

λ, v = eigen(A);

# what are the eigenvalues of `A`?
λ

# The eigenvectors of `A`:
v

# Check
v' * v

#=
# Some equations
=#

#=
Some long equations:
$$
\begin{align*}
A & = \int_0^{\infty} f(x)\,\mathrm{d}x \\
  & = F(x)
\end{align*}
$$
Some texts.
=#

#=
Some code listing
```julia
α = 1.0
β = 2.1
α + β

mutable struct MyStruct
  A::Vector{Matrix{ComplexF64}}
  B::Vector{Matrix{Float64}}
end
```
=#