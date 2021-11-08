using FFTW

function test_inplace!(plan, Ns, A)
    Ar = reshape(A, Ns)
    plan*Ar
    return
end

function G_to_R!( plan, Ns, fG::AbstractArray{ComplexF64} )
    f = reshape(fG, Ns)
    plan*f
    return
end

function mainv1()
    Ns = (50,50,50)
    planfw = plan_fft!(zeros(ComplexF64,Ns))
    A = rand(ComplexF64,prod(Ns))
    println("Before A[1] = ", A[1])
    @time test_inplace!(planfw, Ns, A)
    println("After A[1] = ", A[1])
end
mainv1()

function mainv2()
    Ns = (50,50,50)
    planfw = plan_fft!(zeros(ComplexF64,Ns))
    A = rand(ComplexF64,prod(Ns))
    println("Before A[1] = ", A[1])
    @time G_to_R!(planfw, Ns, A)
    println("After A[1] = ", A[1])
    println("size(A) = ", size(A))
end
mainv2()