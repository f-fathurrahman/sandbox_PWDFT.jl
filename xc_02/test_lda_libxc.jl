using Printf
using PWDFT

include("Libxc_old.jl")

function main()
    func_id = 1 # LDA_X
    #func_id = 7 # LDA_C_VWN

    rho   = [1.00]

    Npoints = length(rho)
    Nspin = 1
    
    # outputs
    exc = zeros(Float64, Npoints)
    vrho = zeros(Float64, Npoints)

    ptr = Libxc_xc_func_alloc()
    Libxc_xc_func_init(ptr, func_id, Nspin)
    
    Libxc_xc_lda_exc_vxc!( ptr, Npoints, rho, exc, vrho )

    Libxc_xc_func_end(ptr)

    for i in 1:Npoints
        @printf("%18.10f: %18.10f %18.10f\n", rho[i], exc[i], vrho[i])
    end

    println("Pass here")

end

main()