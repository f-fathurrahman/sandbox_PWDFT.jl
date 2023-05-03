using Printf
using Random

using PWDFT

Random.seed!(1234)

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")

include(joinpath(DIR_PWDFT, "utilities", "PWSCFInput.jl"))
include(joinpath(DIR_PWDFT, "utilities", "init_Ham_from_pwinput.jl"))

function main(;filename=nothing)

    # Read PWSCF input file and construct a Hamiltonian from it
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    Rhoe, RhoeG = atomic_rho_g(Ham)
    if any(Ham.pspotNL.are_paw)
        PAW_atomic_becsum!(Ham)
    end

    Ehartree, Exc, Evtxc = update_from_rhoe!( Ham, nothing, Rhoe, RhoeG )

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

end

main()
