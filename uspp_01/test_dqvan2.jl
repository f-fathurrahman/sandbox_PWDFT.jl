using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("dqvan2.jl")

function test_main(; filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    atoms = Ham.atoms
    pspotNL = Ham.pspotNL
    pw = Ham.pw

    Ng = pw.gvec.Ng
    G = pw.gvec.G
    G2 = pw.gvec.G2
    dQfuncG = zeros(ComplexF64, Ng)
    nh = pspotNL.nh
    lmaxkb = pspotNL.lmaxkb

    # Input
    isp = 1
    ih = 3
    jh = 3

    @assert ih <= nh[isp]
    @assert jh <= nh[isp]
    @assert isp <= atoms.Nspecies

    # Input
    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    Ylm_real_qe!(_lmax, pw.gvec.G, ylmk0) # Ylm_real_qe accept l value starting from 0

    ipol = 1
    dylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    PWDFT.dYlm_real_qe!(_lmax, pw.gvec.G, dylmk0, ipol)

    dqvan2!(pspotNL, ipol, ih, jh, isp, G, G2, ylmk0, dylmk0, dQfuncG)
    for i in 1:10
        @printf("%5d [%18.10f,%18.10f]\n", i, real(dQfuncG[i]), imag(dQfuncG[i]))
    end

    println("sum dQfuncG = ", sum(abs.(dQfuncG)))
end

test_main()

