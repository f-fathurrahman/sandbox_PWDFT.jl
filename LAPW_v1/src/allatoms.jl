#using ElkDFTWrapper: elk_solve_atom!

function allatoms!( atsp_vars::AtomicSpeciesVars )

    # FIXME
    #xctsp = atsp_vars.xctsp
    xcgrad = false  # hardcoded for now
    xc_calc = LibxcXCCalculator(x_id=1, c_id=12) # for LDA PW92

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

    # Temporary array for storing radial wf
    # (not used outside this function)
    rwf = Vector{Array{Float64,3}}(undef,Nspecies)
    for isp in 1:Nspecies
        rwf[isp] = zeros(Float64,nrsp[isp],2,nstsp[isp])
    end
    #
    for isp in 1:Nspecies
        solve_atom!(
            ptnucl, spzn[isp], nstsp[isp], nsp[isp], lsp[isp], ksp[isp],
            occsp[isp], xc_calc, xcgrad, rsp[isp], evalsp[isp], rhosp[isp],
            vrsp[isp], rwf[isp]
        )
        for ist in 1:nstsp[isp]
            @printf("%3d %18.10f\n", ist, evalsp[isp][ist])
        end
        @printf("sum(rhosp) = %18.10e\n", sum(rhosp[isp]))
    end
    return
end