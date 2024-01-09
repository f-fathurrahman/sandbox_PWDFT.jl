using Printf
using LinearAlgebra
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("gen_neighbor_shells.jl")

function main()
    Ham, pwinput = init_Ham_from_pwinput()
    LatVecs = Ham.pw.LatVecs
    RecVecs = Ham.pw.RecVecs

    dtau = [0.0, 0.0, 0.0]
    mxr = 50
    rmax = 11.752339049876390 # choose the one used in QE test

    r = zeros(Float64, 3, mxr)
    r2 = zeros(Float64, mxr)

    Nvecr = gen_neighbor_shells!( dtau, rmax, LatVecs, RecVecs, r, r2 )    
    println("Nvecr = ", Nvecr)
    for i in 1:Nvecr
        @printf("%5d %18.10f %18.10f %18.10f %18.10f\n", i, r[1,i], r[2,i], r[3,i], sqrt(r2[i]))
    end
    println("rmax = ", rmax)

    return
end

main()
