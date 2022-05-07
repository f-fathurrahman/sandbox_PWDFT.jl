using Printf
using PWDFT: Atoms
# FIXME: Conflict LAPWDFT.XCCalculator and PWDFT.XCCalculator?

using LAPWDFT
using ElkDFTWrapper: elk_solve_atom!


function main(atsymb)

    # Not relevant for atomic problem
    LatVecs = zeros(3,3)
    A = 5.0
    LatVecs[1,:] = [A, 0.0, 0.0]
    LatVecs[2,:] = [0.0, A, 0.0]
    LatVecs[3,:] = [0.0, 0.0, A]

    atoms = Atoms(xyz_string_frac="""
    1

    $atsymb  0.0  0.0  0.0
    """, in_bohr=true, LatVecs=LatVecs)

    Nspecies = atoms.Nspecies
    spsymb = atoms.SpeciesSymbols

    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    for isp in 1:Nspecies
        readspecies!(isp, "DATA_species/"*spsymb[isp]*".in", atsp_vars, mt_vars, apwlo_vars)
    end

    init_zero!(mt_vars)

    checkmt!(atoms, mt_vars )
    genrmesh!(atoms, atsp_vars, mt_vars)
    init_packed_mtr!( mt_vars )


    xcgrad = false  # hardcoded for now
    xc_calc = LibxcXCCalculator()

    spzn = atsp_vars.spzn
    nstsp = atsp_vars.nstsp
    nsp = atsp_vars.nsp
    lsp = atsp_vars.lsp
    ksp = atsp_vars.ksp
    occsp = atsp_vars.occsp
    nrsp = atsp_vars.nrsp
    rsp = atsp_vars.rsp
    evalsp = atsp_vars.evalsp
    ptnucl = atsp_vars.ptnucl
    # Will be calculated in solve_atom!
    rhosp = atsp_vars.rhosp
    vrsp = atsp_vars.vrsp

    Nspecies = length(rsp)

    @assert Nspecies == 1

    # Temporary array for storing radial wf
    # (not used outside this function)
    rwf = Vector{Array{Float64,3}}(undef,Nspecies)
    for isp in 1:Nspecies
        rwf[isp] = zeros(Float64,nrsp[isp],2,nstsp[isp])
    end
    
    # speed of light in atomic units (=1/alpha) (CODATA 2018)
    sol = 137.035999084
    # scaled speed of light
    solsc = sol

    isp = 1

    #println("Calling solve_atom! in allatoms!")
    #println("solsc = ", solsc)
    #println("ptnucl = ", ptnucl)
    #println("spzn = ", spzn[isp])
    #println("nst = ", nstsp[isp])
    #println("n = ", nsp[isp])
    #println("l = ", lsp[isp])
    #println("k = ", ksp[isp])
    #println("occ = ", occsp[isp])
    #println("xc_calc = ", xc_calc)
    #println("xcgrad = ", xcgrad)
    #println("nr = ", nrsp[isp])
    #println("size rsp = ", size(rsp[isp]))
    #println("size evalsp = ", size(evalsp[isp]))
    #println("size rhosp = ", size(rhosp[isp]))
    #println("size vrsp = ", size(vrsp[isp]))
    #println("size rwf = ", size(rwf[isp]))

    solve_atom!(
        solsc, ptnucl, spzn[isp], nstsp[isp], nsp[isp], lsp[isp], ksp[isp],
        occsp[isp], xc_calc, xcgrad, nrsp[isp], rsp[isp], evalsp[isp], rhosp[isp],
        vrsp[isp], rwf[isp]
    )

    #xctype = [3, 0, 0]
    #elk_solve_atom!(
    #    solsc, ptnucl, spzn[isp], nstsp[isp], nsp[isp], lsp[isp], ksp[isp],
    #    occsp[isp], xctype, xcgrad, nrsp[isp], rsp[isp], evalsp[isp], rhosp[isp],
    #    vrsp[isp], rwf[isp]
    #)


    for ist in 1:nstsp[isp]
        @printf("%3d %18.10f\n", ist, evalsp[isp][ist])
    end
    @printf("sum(rhosp) = %18.10e\n", sum(rhosp[isp]))

end

@time main("Pt")
