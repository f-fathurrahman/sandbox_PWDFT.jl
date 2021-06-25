using Printf
using LinearAlgebra
using OffsetArrays
using SpecialFunctions: sphericalbesselj

using PWDFT
using LAPWDFT

import PyPlot
const plt = PyPlot

include("create_atoms.jl")

function create_atsp_mt_apwlo_vars(atoms)

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
    #
    allatoms!(atsp_vars)

    return atsp_vars, mt_vars, apwlo_vars
end


function create_pwgrid(
    atoms, sym_info, mt_vars; rgkmax=7.0, gmaxvr=12.0, kpt_grid=[1,1,1]
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # use average muffin-tin radius (default)
    println("Using average muffin-tin radius to determine gkmax")
    rsum = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        rsum = rsum + mt_vars.rmt[isp]
        println("rmt[isp] = ", mt_vars.rmt[isp])
    end
    rsum = rsum/Natoms
    gkmax = rgkmax/rsum
    println("gkmax = ", gkmax)

    if gmaxvr <= 2*gkmax
        println("gkmax is larger than given/default gmaxvr")
        println("Using gmaxvr = 2*gkmax")
        gmaxvr = 2*gkmax
    end
    println("gmaxvr = ", gmaxvr)

    ecutrho = 0.5*gmaxvr^2
    ecutwfc = 0.5*gkmax^2

    dual = ecutrho/ecutwfc
    println("dual = ", dual)

    pw = PWGrid( ecutwfc, atoms.LatVecs, dual=dual,
        kpoints=KPoints(atoms, kpt_grid, [0,0,0], sym_info.s)
    )
    return pw
end



function main()

    atoms = create_Si_atom()

    atsp_vars, mt_vars, apwlo_vars = create_atsp_mt_apwlo_vars(atoms)

    sym_info = SymmetryInfo(atoms)
    println("sym_info.Nsyms = ", sym_info.Nsyms)
    println("sym_info.Nrots = ", sym_info.Nrots)

    pw = create_pwgrid( atoms, sym_info, mt_vars;
        kpt_grid=[1,1,1], rgkmax=7.0, gmaxvr=12.0
    )
    println(pw)

    # FIXME: need to precompute these?
    sfacg = calc_strfact(atoms, pw)
    ylmg = zeros(ComplexF64, mt_vars.lmmaxo, pw.gvec.Ng)
    genylmg!(mt_vars.lmaxo, pw.gvec.G, ylmg)

    #
    # Initialize rhomt from rhosp
    #

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    nrsp = atsp_vars.nrsp
    rsp = atsp_vars.rsp
    rhosp = atsp_vars.rhosp

    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    ecutrho = pw.ecutrho
    gmaxvr = sqrt(2*ecutrho)
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    idx_g2r = pw.gvec.idx_g2r

    npmt = mt_vars.npmt
    nrmt = mt_vars.nrmt
    npcmt = mt_vars.npcmt
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    rcmt = mt_vars.rcmt
    lmmaxi = mt_vars.lmmaxi
    lmmaxo = mt_vars.lmmaxo

    lmax = min(mt_vars.lmaxi,1)
    epslat = 1e-6

    # Initialize rhomt and rhoir
    Npoints = prod(pw.Ns)
    rhomt = Vector{Vector{Float64}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        rhomt[ia] = zeros(Float64, npmt[isp])
    end

    nrcmtmax = maximum(nrcmt)
    npcmtmax = maximum(npcmt)
    jl = OffsetArray(
        zeros(Float64,lmax+1,nrcmtmax), 0:lmax, 1:nrcmtmax
    )
    zfmt = zeros(ComplexF64,npcmtmax)
    zfft = zeros(ComplexF64,Npoints)
  
    #
    # compute the superposition of all the atomic density tails
    #   
    for isp in 1:Nspecies
        # local arrays inside the loop
        wr = zeros(Float64,nrsp[isp])
        fr = zeros(Float64,nrsp[isp])
        nr = nrmt[isp]
        nrs = nrsp[isp]
        nro = nrs - nr + 1
        # determine the weights for the radial integral
        # integrate for radial points outside the muffin-tin
        idx = nr:nr+nro-1 # or better: idx = nr:nrs
        @views wsplint!( nro, rsp[isp][idx] ,wr[idx] )
        for ig in 1:Ng
            t1 = sqrt(G2[ig])
            # spherical bessel function j_0(x) times the atomic density tail
            if t1 > epslat
                # G != 0 term
                t2 = 1.0/t1
                for ir in nr:nrs
                  x = t1*rsp[isp][ir]
                  fr[ir] = t2*sin(x)*rhosp[isp][ir]*rsp[isp][ir]
                end
            else
                # G=0 term
                for ir in nr:nrs
                    fr[ir] = rhosp[isp][ir] * rsp[isp][ir]^2
                end
            end
            @views t1 = dot( wr[nr:nrs], fr[nr:nrs] )
            # apply low-pass filter
            t1 = t1*exp(-4.0*G2[ig]/gmaxvr^2)
            ffg = (4*pi/CellVolume)*t1
            #
            ip = pw.gvec.idx_g2r[ig]
            zfft[ip] = zfft[ip] + ffg*sfacg[ig,isp]  # do not use conj
        end
    end

    for ia in 1:Natoms
        isp = atm2species[ia]
        nrc = nrcmt[isp]
        nrci = nrcmti[isp]
        irco = nrci + 1
        zfmt[1:npcmt[isp]] .= 0.0
        for ig in 1:Ng
            ip = idx_g2r[ig]
            for irc in 1:nrc
                x = sqrt(G2[ig])*rcmt[isp][irc]
                for l in 0:lmax
                    jl[l,irc] = sphericalbesselj(l, x)
                end
            end
            z1 = 4*pi*zfft[ip] * conj(sfacg[ig,isp]) # XXX using conj
            lm = 0
            for l in 0:lmax
                z2 = im^l * z1
                for m in -l:l
                    lm = lm + 1
                    z3 = z2*conj(ylmg[lm,ig])
                    i = lm
                    for irc in 1:nrci
                        zfmt[i] = zfmt[i] + jl[l,irc]*z3
                        i = i + lmmaxi
                    end
                    for irc in irco:nrc
                        zfmt[i] = zfmt[i] + jl[l,irc]*z3
                        i = i + lmmaxo
                    end
                end
            end
        end
        #println("Before z_to_rf_mt")
        #println(rhomt[ia][1:4])
        z_to_rf_mt!( mt_vars, nrc, nrci, zfmt, rhomt[ia] )
    end

    ia = 1
    isp = atm2species[ia]

    npmti = mt_vars.npmti
    nrmti = mt_vars.nrmti
    idx_outer = npmti[isp]+1:npmt[isp]

    rho_inner = reshape( rhomt[ia][1:npmti[isp]], (lmmaxi,nrmti[isp]) )
    rho_outer = reshape( rhomt[ia][idx_outer], (lmmaxo,nrmt[isp]-nrmti[isp]) )

    plt.clf()
    plt.plot(rho_inner[1,:], marker="o")
    plt.grid(true)
    plt.savefig("IMG_rho_inner_before.pdf")

    plt.clf()
    plt.plot(rho_outer[1,:], marker="o")
    plt.grid(true)
    plt.savefig("IMG_rho_outer_before.pdf")

    # convert the density from a coarse to a fine radial mesh
    rf_mt_c_to_f!( atoms, atsp_vars, mt_vars, rhomt )
    rho_inner = reshape( rhomt[ia][1:npmti[isp]], (lmmaxi,nrmti[isp]) )
    rho_outer = reshape( rhomt[ia][idx_outer], (lmmaxo,nrmt[isp]-nrmti[isp]) )

    plt.clf()
    plt.plot(rho_inner[1,:], marker="o")
    plt.grid(true)
    plt.savefig("IMG_rho_inner_after.pdf")

    plt.clf()
    plt.plot(rho_outer[1,:], marker="o")
    plt.grid(true)
    plt.savefig("IMG_rho_outer_after.pdf")

end

main()