using Printf
using OffsetArrays
using LinearAlgebra

include("../ylm_real/Ylm_real_qe.jl")
#include("../ylm_real/_generate_lm_indices_qe.jl")
include("calc_clebsch_gordan.jl")

function test_main()
    ap, lpx, lpl = calc_clebsch_gordan(1)
end

test_main()

