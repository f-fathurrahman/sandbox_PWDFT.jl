using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

include("gen_us_dj.jl")
include("gen_us_dy.jl")

function test_main(; filename=nothing)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    atoms = Ham.atoms
    pspotNL = Ham.pspotNL
    pw = Ham.pw

    ik = 2
    @assert ik <= Ham.pw.gvecw.kpoints.Nkpt

    Ngwx = maximum(Ham.pw.gvecw.Ngw)
    NbetaNL = Ham.pspotNL.NbetaNL
    dvkb = zeros(ComplexF64, Ngwx, NbetaNL)

    _gen_us_dy!(ik, [1.0, 0.0, 0.0], Ham.atoms, Ham.pw, Ham.pspots, Ham.pspotNL, dvkb)
    println("sum dvkb from _gen_us_dy = ", sum(dvkb))

    _gen_us_dj!(ik, Ham.atoms, Ham.pw, Ham.pspots, Ham.pspotNL, dvkb)
    println("sum dvkb from _gen_us_dj = ", sum(dvkb))

end

test_main()

