mutable struct AtomicSpeciesVars
    # species files path
    sppath::String
    # species filenames
    spfname::Vector{String}
    # species name
    spname::Vector{String}
    # species symbol
    spsymb::Vector{String}
    # species nuclear charge
    spzn::Vector{Float64}
    # ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
    # the nuclei have a finite spherical distribution
    ptnucl::Bool
    # nuclear Coulomb potential
    vcln::Vector{Vector{Float64}}
    # species electronic charge
    spze::Vector{Float64}
    # species mass
    spmass::Vector{Float64}
    # smallest radial point for each species
    rminsp::Vector{Float64}
    # effective infinity for species
    rmaxsp::Vector{Float64}
    # number of radial points to effective infinity for each species
    nrsp::Vector{Int64}
    # maximum allowed states for each species
    maxstsp::Int64 # parameter: 40
    # number of states for each species
    nstsp::Vector{Int64}
    # state principle quantum number for each species
    nsp::Vector{Vector{Int64}} # (maxstsp,maxspecies)
    # state l value for each species
    lsp::Vector{Vector{Int64}} # (maxstsp,maxspecies)
    # state k value for each species
    ksp::Vector{Vector{Int64}} #(maxstsp,maxspecies)
    # spcore is .true. if species state is core
    spcore::Vector{Vector{Bool}} #(maxstsp,maxspecies)
    # total number of core states
    nstcr::Int64
    # state eigenvalue for each species
    evalsp::Vector{Vector{Float64}} # (maxstsp,maxspecies)
    # state occupancy for each species
    occsp::Vector{Vector{Float64}} # (maxstsp,maxspecies)
    # species radial mesh to effective infinity
    rsp::Vector{Vector{Float64}}
    # species charge density
    rhosp::Vector{Vector{Float64}}
    # species self-consistent potential
    vrsp::Vector{Vector{Float64}}
    # exchange-correlation type for atomic species (the converged ground-state of
    # the crystal does not depend on this choice)
    xctsp::Tuple{Int64,Int64,Int64}
end
# rlsp is not used. It is also removed in the future version of Elk
# Some variables are removed: nucl variables and cutoff energies for species generation

function AtomicSpeciesVars( atoms::Atoms, specs_info::Vector{SpeciesInfo} )
    
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Nspecies = length(specs_info)

    sppath = "" # This is not really used for now

    # Copy information from specs_info
    spfname = Vector{String}(undef, Nspecies)
    spname = Vector{String}(undef, Nspecies)
    spsymb = Vector{String}(undef, Nspecies)
    spzn = zeros(Float64, Nspecies)
    spmass = zeros(Float64, Nspecies)
    rminsp = zeros(Float64, Nspecies)
    rmaxsp = zeros(Float64, Nspecies)
    nstsp = zeros(Int64, Nspecies)
    nsp = Vector{Vector{Int64}}(undef,Nspecies)
    lsp = Vector{Vector{Int64}}(undef,Nspecies)
    ksp = Vector{Vector{Int64}}(undef,Nspecies)
    spcore = Vector{Vector{Bool}}(undef,Nspecies)
    evalsp = Vector{Vector{Float64}}(undef,Nspecies)
    occsp  = Vector{Vector{Float64}}(undef,Nspecies)
    for isp in 1:Nspecies
        spfname[isp] = specs_info[isp].filename
        spname[isp] = specs_info[isp].spname
        spsymb[isp] = specs_info[isp].spsymb
        spzn[isp] = specs_info[isp].spzn
        spmass[isp] = specs_info[isp].spmass
        rminsp[isp] = specs_info[isp].rminsp
        rmaxsp[isp] = specs_info[isp].rmaxsp
        nstsp[isp] = specs_info[isp].nstsp
        nsp[isp] = specs_info[isp].nsp
        lsp[isp] = specs_info[isp].lsp
        ksp[isp] = specs_info[isp].ksp
        spcore[isp] = specs_info[isp].spcore
        evalsp[isp] = specs_info[isp].evalsp
        occsp[isp] = specs_info[isp].occsp
    end

    nstcr = 0
    spze = zeros(Float64, Nspecies)
    for ia in 1:Natoms
        isp = atm2species[ia]
        # valence charge
        for ist in 1:nstsp[isp]
            spze[isp] += occsp[isp][ist]
            if spcore[isp][ist] 
                nstcr += 2*ksp[isp][ist]
            end 
        end
    end

    nrsp = zeros(Int64, Nspecies)
    # estimate the number of radial mesh points to infinity
    for isp in 1:Nspecies
        rmt = specs_info[isp].rmt
        nrmt = specs_info[isp].nrmt    
        # logarithmic mesh
        t1 = log(rmaxsp[isp]/rminsp[isp])/log(rmt/rminsp[isp])
        t2 = (nrmt - 1)*t1
        nrsp[isp] = round(Int64,t2) + 1 # XXX compare with nint
    end

    rsp = Vector{Vector{Float64}}(undef, Nspecies)
    rhosp = Vector{Vector{Float64}}(undef, Nspecies)
    vrsp = Vector{Vector{Float64}}(undef, Nspecies)
    for isp in 1:Nspecies
        rsp[isp] = zeros(Float64, nrsp[isp])
        rhosp[isp] = zeros(Float64, nrsp[isp])
        vrsp[isp] = zeros(Float64, nrsp[isp])
    end
    # rhosp and vrsp will be given after in allatoms
    #
    # Setup radial mesh
    for isp in 1:Nspecies
        nrmt = specs_info[isp].nrmt
        rmt = specs_info[isp].rmt
        t1 = 1.0/(nrmt - 1)
        # logarithmic mesh
        t2 = log(rmt/rminsp[isp])
        for ir in 1:nrsp[isp]
            rsp[isp][ir] = rminsp[isp]*exp( (ir-1) * t1 * t2)
        end
    end
    #
    # determine the nuclear Coulomb potential
    ptnucl = true
    # spherical harmonic for l=m=0
    y00 = 0.28209479177387814347
    t1 = 1.0/y00
    #
    vcln = Vector{Vector{Float64}}(undef,Nspecies)
    for isp in 1:Nspecies
        nr = nrsp[isp]
        vcln[isp] = zeros(Float64,nr)
        potnucl!( ptnucl, rsp[isp], spzn[isp], vcln[isp] )
        # scale with 1/y00
        for ir in 1:nr
            vcln[isp][ir] = t1*vcln[isp][ir]
        end
    end

    maxstsp = 40 # parameter, not used?
    xctsp  = (3,0,0) # not used, always defaulting to this?

    return AtomicSpeciesVars(
        sppath, spfname, spname, spsymb, spzn,
        ptnucl, vcln, spze, spmass, rminsp, rmaxsp, nrsp,
        maxstsp, nstsp,
        nsp, lsp, ksp, spcore, nstcr, evalsp, occsp, rsp, rhosp, vrsp, xctsp
    )
end


function AtomicSpeciesVars( Nspecies::Int64 )
    
    sppath = ""
    spfname = Vector{String}(undef, Nspecies)
    spname = Vector{String}(undef, Nspecies)
    spsymb = Vector{String}(undef, Nspecies)
    
    spzn = zeros(Float64, Nspecies)    

    ptnucl = true
    
    vcln = Vector{Vector{Float64}}(undef,Nspecies)
    
    spze = zeros(Float64, Nspecies)
    spmass = zeros(Float64, Nspecies)
    
    rminsp = zeros(Float64, Nspecies)
    rmaxsp = zeros(Float64, Nspecies)
    nrsp = zeros(Int64, Nspecies)
    maxstsp = 40
    nstsp = zeros(Int64, Nspecies)

    nsp = Vector{Vector{Int64}}(undef,Nspecies)
    lsp = Vector{Vector{Int64}}(undef,Nspecies)
    ksp = Vector{Vector{Int64}}(undef,Nspecies)
    spcore = Vector{Vector{Bool}}(undef,Nspecies)

    nstcr = 0
    evalsp = Vector{Vector{Float64}}(undef,Nspecies)
    occsp = Vector{Vector{Float64}}(undef,Nspecies)

    # XXX Will be setup properly later
    rsp = Vector{Vector{Float64}}(undef,Nspecies)
    rhosp = Vector{Vector{Float64}}(undef,Nspecies)
    vrsp = Vector{Vector{Float64}}(undef,Nspecies)

    xctsp  = (3,0,0) # not used, always defaulting to this?

    return AtomicSpeciesVars(
        sppath, spfname, spname, spsymb, spzn,
        ptnucl,
        vcln, spze, spmass, rminsp, rmaxsp, nrsp,
        maxstsp, nstsp,
        nsp, lsp, ksp, spcore, nstcr, evalsp, occsp, rsp, rhosp, vrsp, xctsp
    )
end


# Setting up vcln
function init_nuclear_pot!( atsp_vars::AtomicSpeciesVars )
    # spherical harmonic for l=m=0
    y00 = 0.28209479177387814347

    rsp = atsp_vars.rsp
    nrsp = atsp_vars.nrsp
    Nspecies = size(nrsp,1)
    ptnucl = atsp_vars.ptnucl
    spzn = atsp_vars.spzn

    # determine the nuclear Coulomb potential
    t1 = 1.0/y00
    for isp in 1:Nspecies
        nr = nrsp[isp]
        atsp_vars.vcln[isp] = zeros(Float64,nr)
        potnucl!( ptnucl, rsp[isp], spzn[isp], atsp_vars.vcln[isp] )
        for ir in 1:nr
            atsp_vars.vcln[isp][ir] = t1*atsp_vars.vcln[isp][ir]
        end
    end
    return
end