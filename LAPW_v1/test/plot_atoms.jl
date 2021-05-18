using Printf
using LAPWDFT

import PyPlot
const plt = PyPlot

function create_lattice_vars()
    LatVecs = zeros(3,3)
    A = 10.0
    LatVecs[1,:] = [A, A, 0.0]
    LatVecs[2,:] = [A, 0.0, A]
    LatVecs[3,:] = [0.0, A, A]
    lattice_vars = LatticeVars( LatVecs )
    return lattice_vars
end

function create_atomic_vars(lattice_vars)

    maxatoms = 200
    maxspecies = 8
    atposl = zeros(3,maxatoms,maxspecies)
    
    Nspecies = 1
    Natoms = [1]

    # species 1, atom 1
    atposl[:,1,1] = [0.0, 0.0, 0.0]

    atomic_vars = AtomicVars(Nspecies, Natoms, atposl, lattice_vars)
end

function main()

    latt_vars = create_lattice_vars()
    atm_vars = create_atomic_vars(latt_vars)

    Nspecies = 1
    atsp_vars = AtomicSpeciesVars(Nspecies)
    mt_vars = MuffinTins(Nspecies)
    apwlo_vars = APWLOVars(Nspecies, mt_vars.lmaxapw)

    readspecies!(1, "DATA_species/H.in", atsp_vars, mt_vars, apwlo_vars)

    init_zero!( mt_vars )

    println("before nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    checkmt!( latt_vars, atm_vars, atsp_vars.spsymb, mt_vars )
    genrmesh!( atm_vars, atsp_vars, mt_vars )
    init_packed_mtr!(mt_vars)

    println("after nrsp = ", atsp_vars.nrsp)
    println("atsp_vars.rsp = ", size(atsp_vars.rsp))

    Nspecies = atm_vars.Nspecies

    xctsp = atsp_vars.xctsp
    xcgrad = false

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
    # Will calculated in solve_atom!
    rhosp = atsp_vars.rhosp
    vrsp = atsp_vars.vrsp

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

    for isp in 1:Nspecies
        solve_atom!(
            solsc, ptnucl, spzn[isp], nstsp[isp], nsp[isp], lsp[isp], ksp[isp],
            occsp[isp], xctsp, xcgrad, nrsp[isp], rsp[isp], evalsp[isp], rhosp[isp],
            vrsp[isp], rwf[isp]
        )


        plt.clf()
        for ist in 1:nstsp[isp]

            if atsp_vars.spcore[isp][ist]
                plt.plot(rsp[isp], rwf[isp][:,1,ist], label="st"*string(ist), color="black")
            else
                plt.plot(rsp[isp], rwf[isp][:,1,ist], label="st"*string(ist))
            end

            @printf("%3d (%2d,%2d,%2d) %8.5f %18.10f\n", ist,
                nsp[isp][ist], lsp[isp][ist], ksp[isp][ist], occsp[isp][ist],
                evalsp[isp][ist]
            )
        end
        plt.xlim(atsp_vars.rminsp[isp], mt_vars.rmt[isp])
        plt.grid()
        plt.legend()
        filename = "IMG_atom_"*atsp_vars.spsymb[isp]*"_rwf.pdf"
        plt.savefig(filename)
    end

end

main()
