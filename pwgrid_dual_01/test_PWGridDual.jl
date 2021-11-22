using PWDFT
using Printf
using LinearAlgebra
using FFTW: plan_fft!, plan_ifft!

include("PWGridDual.jl")

function main()
    pw1 = PWGridDual(15.0, gen_lattice_fcc(5.0), dual=5.0)
    println(pw1)
    println("Size (KB) = ", Base.summarysize(pw1)/1024)

    pw2 = PWGridDual(15.0, gen_lattice_fcc(5.0))
    println(pw2)
    println("Size (KB) = ", Base.summarysize(pw2)/1024)
end

main()