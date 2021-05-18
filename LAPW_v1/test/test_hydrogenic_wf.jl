push!(LOAD_PATH, pwd())

using Printf
using LinearAlgebra
using PWDFT
using LAPWDFT

import PyPlot
const plt = PyPlot

function main()

    LatVecs = zeros(3,3)
    A = 10.0
    LatVecs[1,:] = [A, 0.0, 0.0]
    LatVecs[2,:] = [0.0, A, 0.0]
    LatVecs[3,:] = [0.0, 0.0, A]

    atoms = Atoms(xyz_string="""
    1

    H  0.0  0.0  0.0
    """, in_bohr=true, LatVecs=LatVecs)

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

    r = atsp_vars.rsp[1]
    Nr = atsp_vars.nrsp[1]
    println("Nr = ", Nr)
    n = 1
    l = 0
    Z = 1.0
    Zr = Z*r
    fr_1s = (Z)^1.5 .* 2 .* exp.(-Zr)
    fr_2s = (Z/2)^1.5 .* (2 .- Zr) .* exp.(-0.5*Zr)
    fr_3s = (Z/3)^1.5 .* 2 .* (1 .- (2/3).*Zr .+ 2/27 .* (Zr).^2) .* exp.(-Zr/3)
    #
    wr = zeros(Float64,Nr)
    wsplint!( Nr, r, wr )
    #
    ss = dot(r.^2 .* fr_1s.^2, wr)
    println("Norm 1s = ", ss)
    #
    ss = dot(r.^2 .* fr_2s.^2, wr)
    println("Norm 2s = ", ss)
    #
    ss = dot(r.^2 .* fr_3s.^2, wr)
    println("Norm 3s = ", ss)

    plt.clf()
    plt.plot(r, fr_1s, label="1s")
    plt.plot(r, fr_2s, label="2s")
    plt.plot(r, fr_3s, label="3s")
    plt.grid(true)
    plt.xlim(0.0, 20.0)
    plt.savefig("IMG_H_wfcs.pdf")
end

main()
