using Printf
using LinearAlgebra
using Random
using PWDFT: ortho_sqrt!

function main()

    Random.seed!(1234)

    Nstates = 4
    Nbasis = 10
    O = zeros(ComplexF64,Nstates,Nstates)
    U = zeros(ComplexF64,Nstates,Nstates)
    C = zeros(ComplexF64,Nstates,Nstates)

    psis_SC = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(psis_SC)
    phi_m0  = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(phi_m0)
    phi_m1  = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(phi_m1)
    phi_m2  = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(phi_m2)
    phi_m3  = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(phi_m3)
    phi_m4  = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(phi_m4)
    phi_m5  = rand(ComplexF64,Nbasis,Nstates); ortho_sqrt!(phi_m5)

    psis = zeros(ComplexF64,Nbasis,Nstates)

    dt = 1.0
    # XL-BOMD parameters
    κ = 1.82
    ω2 = κ/dt^2
    α = 0.018
    c0 = -6.0
    c1 = 14.0
    c2 = -8.0
    c3 = -3.0
    c4 = 4.0
    c5 = -1.0

    # Alignment w.r.t phi_m5 (the oldest wavefunction)
    do_align!(phi_m4, phi_m5)
    do_align!(phi_m3, phi_m5)
    do_align!(phi_m2, phi_m5)
    do_align!(phi_m1, phi_m5)
    do_align!(phi_m0, phi_m5)

    display_orthonormality("phi_m5", phi_m5)
    display_orthonormality("phi_m4", phi_m4)
    display_orthonormality("phi_m3", phi_m3)
    display_orthonormality("phi_m2", phi_m2)
    display_orthonormality("phi_m1", phi_m1)
    display_orthonormality("phi_m0", phi_m0)

    @views O[:,:] = psis_SC' * phi_m5
    @views U[:,:] = inv(sqrt(O*O')) * O

    println("O before do_align! = ")
    display(O); println()
    println("U before do_align! = ")
    display(U); println()
    
    do_align!(psis_SC, phi_m5)

    @views O[:,:] = psis_SC' * phi_m5
    @views U[:,:] = inv(sqrt(O*O')) * O

    println("O after do_align! = ")
    display(O); println()
    println("U after do_align! = ")
    display(U); println()

    # Alignment
    @views O[:,:] = psis_SC' * phi_m0
    @views U[:,:] = inv(sqrt(O*O')) * O
    @views psis[:,:] = 2*phi_m0 - phi_m1 +
        κ * ( psis_SC*U - phi_m0 ) +
        α * (
            c0*phi_m0 + c1*phi_m1 +
            c2*phi_m2 + c3*phi_m3 +
            c4*phi_m4 + c5*phi_m5
        )

    display_orthonormality("psis", psis)

end

function do_align!(psi, psi0)
    O = psi' * psi0
    U = inv(sqrt(O*O')) * O
    psi[:,:] = psi[:,:]*U
    return
end

function display_orthonormality(s::String, psi::Matrix{ComplexF64}; ist_ref=2)
    Nstates = size(psi,2)
    println(s)
    for ist in 1:Nstates
        s = dot(psi[:,ist], psi[:,ist_ref])
        @printf("dot(%3d,%3d) = %18.10f + %18.10f\n", ist, ist_ref, real(s), imag(s))
    end
    return
end

main()