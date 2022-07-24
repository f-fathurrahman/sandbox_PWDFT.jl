using Printf
using OffsetArrays
using SpecialFunctions: sphericalbesselj
using Random

using PWDFT

include("../pwscf_02/PWSCFInput.jl")
include("../pwscf_02/init_Ham_from_pwinput.jl")

include("atomic_rho_g.jl")
include("dense_to_smooth.jl")
include("update_from_rhoe.jl")
include("newd.jl")
include("op_S.jl")
include("calc_rhoe_uspp.jl")

#atoms, pw, pspots = init_test_main()
function test_main()

    Random.seed!(1234)

    Ham = init_Ham_from_pwinput()

    Rhoe, RhoeG = atomic_rho_g(Ham)

    update_from_rhoe!( Ham, Rhoe, RhoeG )

    calc_newDeeq!( Ham ) # should be included in update_from_rhoe!

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    
    Nspin = Ham.electrons.Nspin
    @assert Nspin == 1

    psiks = Vector{Matrix{ComplexF64}}(undef,Nkpt)
    for ik in 1:Nkpt
        Ham.ik = ik # set the current k index of H
        Ngwk = Ham.pw.gvecw.Ngw[ik]
        psiks[ik] = randn(ComplexF64, Ngwk, Nstates)
        O = psiks[ik]' * op_S(Ham, psiks[ik])
        U = inv(sqrt(O))
        psiks[ik][:,:] = psiks[ik] * U
    end

    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints

    Rhoe = zeros(Float64, Npoints, Nspin)
    calc_rhoe_uspp!(Ham, psiks, Rhoe )
    println("integ Rhoe = ", sum(Rhoe)*dVol)

    nhm = Ham.pspotNL.nhm
    Natoms = Ham.atoms.Natoms
    ik = 1
    ispin = 1
    Nbecsum = Int64( nhm * (nhm + 1)/2 )
    becsum = zeros(Float64, Nbecsum, Natoms, Nspin)

    _add_becsum!(ik, ispin, Ham, psiks, becsum)

    _add_usdens!(Ham, becsum, Rhoe)

    println("integ Rhoe = ", sum(Rhoe)*dVol)

end

test_main()
