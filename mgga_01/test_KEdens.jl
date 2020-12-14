using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("../get_default_psp.jl")

function do_calc(molname)

    Random.seed!(1234)

    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", molname*".xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms, xcfunc="PBE")
    
    ecutwfc = 15.0
    
    Ham = Hamiltonian(atoms, pspfiles, ecutwfc, xcfunc="SCAN")
    psiks = rand_BlochWavefunc(Ham)

    Rhoe = calc_rhoe(Ham, psiks)

    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints
    KEdens = zeros(Float64, Npoints)

    calc_KEdens!(1, Ham.pw, psiks[1], KEdens)
    println("integ KEdens = ", sum(KEdens)*dVol)

    epsxc = calc_epsxc_SCAN( Ham.xc_calc, Ham.pw, psiks, Rhoe[:,1] )
    println("E_xc = ", dot(epsxc, Rhoe)*dVol)

    V_xc = zeros(Float64,Npoints)
    calc_Vxc_SCAN!( Ham.xc_calc, Ham.pw, psiks, Rhoe[:,1], V_xc )
    println("integ V_xc = ", sum(V_xc)*dVol)
    println("sub abs V_xc = ", sum(abs.(V_xc)))
    println("sub abs Vtau = ", sum(abs.(Ham.xc_calc.Vtau)))

    Vpsiks = zeros_BlochWavefunc(Ham)
    op_Vtau!(Ham, psiks, Vpsiks)
    println("dot psiks, Vpsiks = ", dot(psiks,Vpsiks))
    println("dot Vpsiks, Vpsiks = ", dot(Vpsiks,Vpsiks))

    println("Pass here")
end


function main()
    Nargs = length(ARGS)
    if Nargs >= 1
        molname = ARGS[1]
    else
        molname = "H2O"
    end
    do_calc(molname)
end

main()
