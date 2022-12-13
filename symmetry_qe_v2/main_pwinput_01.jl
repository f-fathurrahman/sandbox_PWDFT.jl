using Printf
using OffsetArrays
using LinearAlgebra
import Serialization
import Random

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

function print_matrix(A, Nrows, Ncols)
    for i in 1:Nrows
        for j in 1:Ncols
            @printf("%10.5f", A[i,j])
        end
        println()
    end
    return
end

function main()
    Ham = init_Ham_from_pwinput()

    sym_info = Ham.sym_info
    Nsyms = sym_info.Nsyms
    sr = sym_info.sr

    println("Nsyms = ", Nsyms)
    for isym in 1:Nsyms
        println("\nisym = ", isym)
        #display(sr[:,:,isym]); println()
        print_matrix(sr[:,:,isym], 3, 3)
    end
end

main()
