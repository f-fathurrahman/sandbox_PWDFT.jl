using Printf
using PWDFT

function _calc_epsxc_Vxc_PBE!(
    Rhoe::AbstractVector{Float64},
    gRhoe2::AbstractVector{Float64},
    epsxc::AbstractVector{Float64},
    V_xc::AbstractVector{Float64},
    Vg_xc::AbstractVector{Float64}
)
    Npoints = size(Rhoe, 1)
    Nspin = 1

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)
 
    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)
 
    ptr = PWDFT.Libxc_xc_func_alloc()
    # exchange part
    PWDFT.Libxc_xc_func_init(ptr, 101, Nspin)
    PWDFT.Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    PWDFT.Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, eps_x, V_x, Vg_x)
    PWDFT.Libxc_xc_func_end(ptr)

    #
    # correlation part
    PWDFT.Libxc_xc_func_init(ptr, 130, Nspin)
    PWDFT.Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    PWDFT.Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, eps_c, V_c, Vg_c)
    PWDFT.Libxc_xc_func_end(ptr)

    #
    PWDFT.Libxc_xc_func_free(ptr)

    # update the outputs:
    @views epsxc[:] .= eps_x[:] .+ eps_c[:]
    @views V_xc[:] .= V_x[:] .+ V_c[:]
    @views Vg_xc[:] = Vg_x[:] + Vg_c[:]  # for gradient correction

    # V_xc need to be corrected later using Vg_xc
    return
end


function _driver_gga!(arho, grhoe2, sxc, v1xc, v2xc)
    Nrmesh = size(arho, 1)
    for ir in 1:Nrmesh
        eex, vvx = PWDFT.XC_x_slater( arho[ir,1] )
        eec, vvc = PWDFT.XC_c_pw( arho[ir,1] )
        #
        sx, v1x, v2x = PWDFT.XC_x_pbe( arho[ir,1], grhoe2[ir] )
        sc, v1c, v2c = PWDFT.XC_c_pbe( arho[ir,1], grhoe2[ir] )
        #
        sxc[ir] = sx + sc + (eex + eec)*arho[ir,1]
        v1xc[ir,1] = v1x + v1c + (vvx + vvc)
        v2xc[ir,1] = v2x + v2c
    end
    return
end



function test_main( arho::Float64, gradx::Vector{Float64} )
    #
    grhoe2 = gradx[1]^2 + gradx[2]^2 + gradx[3]^2

    # Array version
    Rhoe = [arho]
    gRhoe2 = [grhoe2]

    epsxc = [0.0]
    V_xc = [0.0]
    Vg_xc =[0.0]
    _driver_gga!(Rhoe, gRhoe2, epsxc, V_xc, Vg_xc)
    @printf("\nInternal version\n")
    @printf("epsxc = %18.10f\n", epsxc[1])
    @printf("V_xc  = %18.10f\n", V_xc[1])
    @printf("Vg_xc = %18.10f\n", Vg_xc[1])

    epsxc = [0.0]
    V_xc = [0.0]
    Vg_xc =[0.0]
    _calc_epsxc_Vxc_PBE!(Rhoe, gRhoe2, epsxc, V_xc, Vg_xc)
    @printf("\nLibxc version\n")
    @printf("epsxc = %18.10f\n", epsxc[1]*arho[1])
    @printf("V_xc  = %18.10f\n", V_xc[1])
    @printf("Vg_xc = %18.10f\n", Vg_xc[1]*2)
    # Please note some conversion needed for epsxc and Vg_xc
    # The correction also can be applied to internal version if
    # Libxc convention is used as reference

end

#test_main(5.0, [2.0, 3.0, 4.0])
#test_main(1.0, [2.0, 3.0, 4.0])
#test_main(1e-10, [2.0, 3.0, 4.0])
#test_main(5.0, 1e-10*[1, 1, 1])
#test_main(5.0, 1e-6*[1, 1, 1])
#test_main(5.0, 1e-2*[1, 1, 1])
#test_main(5.0, 0.1*[1, 1, 1])
test_main(1.0, 0.1*[1, 1, 1])