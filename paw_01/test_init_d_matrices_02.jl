using Printf
using OffsetArrays
using LinearAlgebra
import Serialization
import Random

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

include("init_d_matrices.jl")
include("print_matrix.jl")

function main()
 
    Ham = init_Ham_from_pwinput()
    
    sym_info = Ham.sym_info
    Nsyms = sym_info.Nsyms
    # These matrices are initialized in pwscf if we are using PAW
    # However, we don't need the PAW info for this
    dy1, dy2, dy3 = init_d_matrices( sym_info.sr )

    println("Nsyms = ", sym_info.Nsyms)
    for isym in 1:Nsyms
        println("\nisym = ", isym)

        println("\ndy1 = ")
        print_matrix(dy1[:,:,isym], 3, 3)

        println("\ndy2 = ")
        print_matrix(dy2[:,:,isym], 5, 5)

        println("\ndy3 = ")
        print_matrix(dy3[:,:,isym], 7, 7)
    end
    
end

main()
