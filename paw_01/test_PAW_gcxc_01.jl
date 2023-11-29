using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

#=
function calc_epsxc_Vxc_PBE!(
    Rhoe::AbstractVector{Float64},
    gRhoe2::AbstractVector{Float64},
    epsxc::AbstractVector{Float64},
    Vxc::AbstractVector{Float64}
)

    xc_calc = Ham.xc_calc

    Npoints = size(Rhoe, 1)
    Nspin = 1

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)
 
    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)
 
    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 101, Nspin)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, eps_x, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 130, Nspin)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, eps_c, V_c, Vg_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    @views epsxc[:] .= eps_x[:] .+ eps_c[:]

    @views Vxc[:] .= V_x[:] .+ V_c[:] # update V_xc (the output)

    # gradient correction
    Vg_xc = Vg_x + Vg_c
    hx = zeros(ComplexF64, pw.Ns)
    hy = zeros(ComplexF64, pw.Ns)
    hz = zeros(ComplexF64, pw.Ns)
    for ip in 1:Npoints # using linear indexing
        hx[ip] = Vg_xc[ip] * gRhoe[1,ip]
        hy[ip] = Vg_xc[ip] * gRhoe[2,ip]
        hz[ip] = Vg_xc[ip] * gRhoe[3,ip]
    end
    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, hx, hy, hz )
    #
    for ip in 1:Npoints
        Vxc[ip] = Vxc[ip] - 2.0*divh[ip]
    end

    return
end
=#

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

function main(;filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    becsum = PAW_atomic_becsum(Ham.atoms, Ham.pspots, Ham.pspotNL, Nspin=1)
    
    println("sum becsum before PAW_symmetrize: ", sum(becsum))
    PAW_symmetrize!(Ham, becsum)
    println("sum becsum after PAW_symmetrize: ", sum(becsum))

    Nspin = 1
    ia = 1 # atom index
    isp = Ham.atoms.atm2species[ia] # species index
    Nrmesh = Ham.pspots[isp].Nr
    l2 = (Ham.pspots[isp].lmax_rho + 1)^2

    # Calculate rho_lm
    rho_lm = zeros(Float64, Nrmesh, l2, Nspin)
    AE = true
    PAW_rho_lm!(AE, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm)
    println("sum rho_lm = ", sum(rho_lm))

    rho_rad = zeros(Float64, Nrmesh, Nspin)
    grad = zeros(Float64, Nrmesh, 3, Nspin) # gradient (r, ϕ, θ)
    grad2 = zeros(Float64, Nrmesh, Nspin) # square modulus of gradient

    atoms = Ham.atoms
    pspots = Ham.pspots
    pspotNL = Ham.pspotNL
    nx = pspotNL.paw.spheres[isp].nx

    if AE
        rho_core = pspots[isp].paw_data.ae_rho_atc
    else
        rho_core = pspots[isp].rho_atc
    end

    r = Ham.pspots[isp].r
    r2 = r.^2

    arho = zeros(Float64, Nrmesh) # 2nd dim removed 
    grhoe2 = zeros(Float64, Nrmesh)


    # These arrays should depend on spin
    v1xc = zeros(Float64, Nrmesh, Nspin)
    v2xc = zeros(Float64, Nrmesh, Nspin)
    # These arrays don't depend on spin
    sxc = zeros(Float64, Nrmesh)

    # for energy, it will be summed over for all nx
    e_rad = zeros(Float64, Nrmesh)
    # for potential, will be processed later, depend on nx
    gc_rad = zeros(Float64, Nrmesh, nx, Nspin)
    h_rad = zeros(Float64, Nrmesh, 3, nx, Nspin)

    energy = 0.0 # This should be accumulated for all ix
    spheres = pspotNL.paw.spheres

    #ix = 1
    #@assert ix <= nx

    for ix in 1:nx

        PAW_lm2rad!(ia, ix, atoms, pspots, pspotNL, rho_lm, rho_rad)
        PAW_gradient!(ia, ix, atoms, pspots, pspotNL,
            rho_lm, rho_rad, rho_core,
            grad2, grad
        )

        for ir in 1:Nrmesh
            arho[ir,1] = abs(rho_rad[ir,1]/r2[ir] + rho_core[ir])
            grhoe2[ir] = grad[ir,1]^2 + grad[ir,2]^2 + grad[ir,3]^2
        end

        _driver_gga!(arho, grhoe2, sxc, v1xc, v2xc)

        # radial stuffs
        for ir in 1:Nrmesh
            e_rad[ir] = sxc[ir] * r2[ir]
            gc_rad[ir,ix,1]  = v1xc[ir,1]
            @views h_rad[ir,1:3,ix,1] = v2xc[ir,1]*grad[ir,1:3,1]*r2[ir]
        end
    
        # integrate to obtain the energy
        energy += PWDFT.integ_simpson(Nrmesh, e_rad, pspots[isp].rab)*spheres[isp].ww[ix]

    end

    println("energy for all ix = ", energy)



    lmax_loc = Ham.pspots[isp].lmax_rho + 1
    gc_lm = zeros(Float64, Nrmesh, l2, Nspin)
    # convert the first part of the GC correction back to spherical harmonics
    PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc, gc_rad, gc_lm)

    # trick to get faster convergence w.r.t to θ
    for ix in 1:nx
        @views h_rad[:,3,ix,1] = h_rad[:,3,ix,1]/ spheres[isp].sin_th[ix]
    end

    # We need the gradient of H to calculate the last part of the exchange
    # and correlation potential. First we have to convert H to its Y_lm expansion
    #PAW_rad2lm!( i, h_rad, h_lm, i%l+rad(i%t)%ladd, nspin_gga )
    lmax_loc_add = lmax_loc + spheres[isp].ladd
    h_lm = zeros(Float64, Nrmesh, 3, lmax_loc_add^2, Nspin)
    
    @views PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc_add, h_rad[:,1,:,:], h_lm[:,1,:,:])
    @views PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc_add, h_rad[:,2,:,:], h_lm[:,2,:,:])
    @views PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc_add, h_rad[:,3,:,:], h_lm[:,3,:,:])

    #PAW_rad2lm3!(ia, atoms, pspotNL, lmax_loc_add, h_rad, h_lm)

    println("lmax_loc = ", lmax_loc)
    println("ladd = ", spheres[isp].ladd)
    println("lmax_loc_add = ", lmax_loc_add)
    
    println("sum abs h_rad = ", sum(abs.(h_rad)))
    println("sum abs h_lm = ", sum(abs.(h_lm)))

    println("h_lm 1 = ", h_lm[1,1,1,1])
    println("h_lm 2 = ", h_lm[1,2,1,1])
    println("h_lm 3 = ", h_lm[1,3,1,1])

    println("max abs h_rad 1 = ", maximum(abs.(h_rad[:,1,:,:])))
    println("max abs h_rad 2 = ", maximum(abs.(h_rad[:,2,:,:])))
    println("max abs h_rad 3 = ", maximum(abs.(h_rad[:,3,:,:])))
    println()
    println("max abs h_lm 1 = ", maximum(abs.(h_lm[:,1,:,:])))
    println("max abs h_lm 2 = ", maximum(abs.(h_lm[:,2,:,:])))
    println("max abs h_lm 3 = ", maximum(abs.(h_lm[:,3,:,:])))

    println("sum sin_th = ", sum(spheres[isp].sin_th))

    div_h = zeros(Float64, Nrmesh, lmax_loc^2, Nspin)
    PAW_divergence!(
        ia, atoms, pspots, pspotNL,
        h_lm, div_h, lmax_loc_add, lmax_loc
    )
    println("sum div_h = ", sum(div_h))

    vout_lm = zeros(Float64, Nrmesh, l2, Nspin)
    # Finally sum it back into v_xc
    for ispin in 1:Nspin
        for lm in 1:l2
            @views vout_lm[1:Nrmesh,lm,ispin] .+= gc_lm[1:Nrmesh,lm,ispin] .- div_h[1:Nrmesh,lm,ispin]
        end
    end
    println("sum abs vout_lm = ", sum(abs.(vout_lm)))

end

main()