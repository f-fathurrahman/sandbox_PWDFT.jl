using Printf
using PWDFT

include("Libxc_old.jl")

function main()
    func_id = 263 # MGGA_X_SCAN
    #func_id = 267 # MGGA_C_SCAN

    #rho   = [1.01, 1.02, 0.03, 1.04, 1.05]
    #sigma = [1.01, 0.02, 0.03, 0.04, 0.05]
    #lapl  = [2.1, 2.2, 2.3, 2.4, 2.6]
    #tau   = [1.1, 2.1, 3.1, 4.1, 5.1]

    rho   = [1.1]
    sigma = [0.0]
    lapl  = [0.0]
    tau   = [0.0]

    Npoints = length(rho)
    Nspin = 1
    # outputs
    exc = zeros(Float64, Npoints)
    vrho = zeros(Float64, Npoints)
    vsigma = zeros(Float64, Npoints)
    vlapl = zeros(Float64, Npoints)
    vtau = zeros(Float64, Npoints)

    ptr = Libxc_xc_func_alloc()
    Libxc_xc_func_init(ptr, func_id, Nspin)
    
    Libxc_xc_mgga_exc_vxc!( ptr, Npoints, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau )

    Libxc_xc_func_end(ptr)

    for i in 1:Npoints
        @printf("%18.10f: %18.10f %18.10f %18.10f %18.10f %18.10f\n",
            rho[i], exc[i], vrho[i], vsigma[i], vlapl[i], vtau[i])
    end

    println("Pass here")

end

main()