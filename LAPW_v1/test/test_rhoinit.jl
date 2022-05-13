using Printf
using LinearAlgebra
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT: Atoms, SymmetryInfo, KPoints, PWGrid
using PWDFT: calc_strfact, G_to_R!, R_to_G!
using LAPWDFT

import PyPlot
const plt = PyPlot

include("create_atoms.jl")
include("debug_rhoinit.jl")

function main()

    #atoms = create_H2O()
    atoms = create_Si_fcc()
    #atoms = create_SiPt_fcc()
    #atoms = create_Si_atom()
    #atoms = create_Pt_atom()

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    for isp in 1:Nspecies
        readspecies!(isp, "DATA_species/"*spsymb[isp]*".in", atsp_vars, mt_vars, apwlo_vars)
    end

    init_zero!( mt_vars )
    checkmt!( atoms, mt_vars )
    genrmesh!( atoms, atsp_vars, mt_vars )
    init_packed_mtr!( mt_vars )

    allatoms!(atsp_vars)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Default value for rgkmax and gmaxvr
    rgkmax = 7.0
    gmaxvr = 12.0

    # use average muffin-tin radius (default)
    rsum = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        rsum = rsum + mt_vars.rmt[isp]
        println("rmt[isp] = ", mt_vars.rmt[isp])
    end
    rsum = rsum/Natoms
    gkmax = rgkmax/rsum

    println("rsum = ", rsum)

    if gmaxvr <= 2*gkmax
        println("Using gmaxvr = 2*gkmax")
        gmaxvr = 2*gkmax
    end
    println("gmaxvr = ", gmaxvr)
    println("gkmax = ", gkmax)

    ecutrho = 0.5*gmaxvr^2
    ecutwfc = 0.5*gkmax^2

    dual = ecutrho/ecutwfc
    println("dual = ", dual)

    sym_info = SymmetryInfo(atoms)
    println("sym_info.Nsyms = ", sym_info.Nsyms)
    println("sym_info.Nrots = ", sym_info.Nrots)

    pw = PWGrid( ecutwfc, atoms.LatVecs, dual=dual,
        kpoints=KPoints(atoms, [2,2,2], [0,0,0], sym_info.s)
    )
    println(pw)


    # Initialize rhomt and rhoir
    npmt = mt_vars.npmt
    Npoints = prod(pw.Ns)
    rhomt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        rhomt[ia] = zeros(Float64, npmt[isp])
    end
    rhoir = zeros(Float64,Npoints)
    #
    #debug_rhoinit!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir )
    rhoinit!( atoms, atsp_vars, mt_vars, pw, rhomt, rhoir )
end

main()
