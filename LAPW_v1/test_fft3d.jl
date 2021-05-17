using PWDFT

function main()

    pw = PWGrid( 10.0, gen_lattice_sc(16.0) )
    Ns = pw.Ns
    Npoints = prod(Ns)

    println("Ns = ", Ns)

    zfft = zeros(ComplexF64, Npoints)
    zfft[:] .= 1.0 + im*2.1
    zfft[1] = 99.0 + im*99.0
    
    println("sum zfft before = ", sum(zfft))
    
    G_to_R!(pw, zfft)

    println("sum zfft after = ", Npoints*sum(zfft))
end

main()