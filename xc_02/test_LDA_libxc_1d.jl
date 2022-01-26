using Printf

include("XCCalculator.jl")
include("Libxc_old.jl")

function main()
  
    Npoints = 1
    Nspin = 1
    Rhoe = [1.1]
    eps_c = zeros(Float64, Npoints)

    ptr = Libxc_xc_func_alloc()
    Libxc_xc_func_init(ptr, 18, Nspin) # LDA_C_1D_CSC
    Libxc_xc_lda_exc!(ptr, Npoints, Rhoe, eps_c)
    Libxc_xc_func_end(ptr)
    Libxc_xc_func_free(ptr)

    @printf("eps_c = %18.10f\n", eps_c[1])
end

main()

