using Printf
using OffsetArrays
using LinearAlgebra
using Serialization: serialize

using Random
Random.seed!(1234)

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("radial_gradient_coarse.jl")
include("spline_stuffs.jl")

function main(;filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Nspecies = Ham.atoms.Nspecies
    for isp in 1:Nspecies
        println(Ham.pspots[isp])
    end

    becsum = PAW_atomic_becsum(Ham.atoms, Ham.pspots, Ham.pspotNL, Nspin=1)
    #becsum[:] .= 1.0
    println("sum becsum before PAW_symmetrize: ", sum(becsum))
    PAW_symmetrize!(Ham, becsum)
    println("sum becsum after PAW_symmetrize: ", sum(becsum))

    withcore = false
    Nspin = 1
    atoms = Ham.atoms
    pspotNL = Ham.pspotNL
    Natoms = atoms.Natoms

    grid_R = PWDFT.init_grid_R( Ham.pw.Ns, Ham.pw.LatVecs )
    posi = zeros(Float64, 3)

    ispin = 1 # spin index

    rho_R = zeros(Float64, size(grid_R,2))

    for ia in 1:Natoms

        isp = atoms.atm2species[ia]
        Nrmesh = Ham.pspots[isp].Nr
        lmax_rho = Ham.pspots[isp].lmax_rho
        lmax_loc = lmax_rho + 1
        l2 = lmax_loc^2
        rho_core = Ham.pspots[isp].paw_data.ae_rho_atc
        r = Ham.pspots[isp].r
        r2 = r.^2

        rho_lm = zeros(Float64, Nrmesh, l2, Nspin)
        rho_lm_ae = zeros(Float64, Nrmesh, l2, Nspin)
        rho_lm_ps = zeros(Float64, Nrmesh, l2, Nspin)
        wsp_lm = zeros(Float64, Nrmesh, l2, Nspin) # for spline

        #
        PAW_rho_lm!(true, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm_ae)
        PAW_rho_lm!(false, ia, Ham.atoms, Ham.pspots, Ham.pspotNL, becsum, rho_lm_ps)
        for ispin in 1:Nspin
            for lm in 1:l2, ir in 1:Nrmesh
                rho_lm[ir,lm,ispin] = ( rho_lm_ae[ir,lm,ispin] - rho_lm_ps[ir,lm,ispin] ) / r2[ir]
            end
            # add core charge to (lm=0 component)
            if withcore
                for ir in 1:Nrmesh
                    rho_lm[ir,1,is] += rho_core[ir] / Nspin * (2.0 * sqrt(pi))
                end
            end
        end

        println("sum rho_lm = ", sum(rho_lm))

        kkbeta = Ham.pspots[isp].kkbeta
        d1y = zeros(Float64, kkbeta)
        d2y = zeros(Float64, kkbeta)

        #println("kkbeta = ", kkbeta)
        #println("Nrmesh = ", Nrmesh)
  
        for lm in 1:l2
            @views radial_gradient_coarse!( r[1:kkbeta], rho_lm[1:kkbeta,lm,ispin], d1y )
            @views radial_gradient_coarse!( r[1:kkbeta], d1y, d2y )
            first  = d1y[1] # first derivative in first point
            second = d2y[1] # second derivative in first point
            #println("first = ", first)
            #println("second = ", second)
            @views init_spline!( r, rho_lm[:,lm,ispin], first, second, wsp_lm[:,lm,ispin] )
            println("sum wsp_lm[:,lm,ispin] = ", sum(wsp_lm[:,lm,ispin]))
        end

        ip = 1

        @views posi[:] = grid_R[:,ip] - atoms.positions[:,ia]
        # FIXME: need to account for periodicity
        distsq = sum(posi.^2)
        println("distsq = ", distsq)

        # prepare spherical harmonics
        ylm_posi = zeros(Float64, l2)
        Ylm_real_qe!(lmax_rho, posi, ylm_posi)

        # Loop over spin
        rho_of_r = 0.0
        for lm in 1:l2
            # do interpolation - distsq depends upon ir
            @views rho = spline_interp( r , rho_lm[:,lm,ispin], wsp_lm[:,lm,ispin], sqrt(distsq) )
            println("rho = ", rho)
            rho_of_r += ylm_posi[lm] * rho        
        end
        println("rho_of_r = ", rho_of_r)
        rho_R[ip] += rho_of_r
    end

    println("sum rho_R = ", sum(rho_R))

end


main()