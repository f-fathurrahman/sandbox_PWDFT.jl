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
    # nuclear radius
    rnucl::Vector{Float64}
    # nuclear volume
    volnucl::Vector{Float64}
    # number of radial mesh points to nuclear radius
    nrnucl::Vector{Int64}
    # number of coarse radial mesh points to nuclear radius
    nrcnucl::Vector{Int64}
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
    # core-valence cut-off energy for species file generation
    ecvcut::Float64
    # semi-core-valence cut-off energy for species file generation
    esccut::Float64
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

function AtomicSpeciesVars( Nspecies::Int64 )
    
    sppath = ""
    spfname = Vector{String}(undef, Nspecies)
    spname = Vector{String}(undef, Nspecies)
    spsymb = Vector{String}(undef, Nspecies)
    
    spzn = zeros(Float64, Nspecies)    

    ptnucl = true

    rnucl = zeros(Float64, Nspecies)

    volnucl = zeros(Float64, Nspecies)

    nrnucl = zeros(Int64, Nspecies)
    nrcnucl = zeros(Int64, Nspecies)
    
    vcln = Vector{Vector{Float64}}(undef,Nspecies)
    
    spze = zeros(Float64, Nspecies)
    spmass = zeros(Float64, Nspecies)
    
    rminsp = zeros(Float64, Nspecies)
    rmaxsp = zeros(Float64, Nspecies)
    nrsp = zeros(Int64, Nspecies)
    maxstsp = 40
    nstsp = zeros(Int64, Nspecies)
    
    ecvcut = 0.0
    esccut = 0.0

    nsp = Vector{Vector{Int64}}(undef,Nspecies)
    lsp = Vector{Vector{Int64}}(undef,Nspecies)
    ksp = Vector{Vector{Int64}}(undef,Nspecies)
    spcore = Vector{Vector{Bool}}(undef,Nspecies)

    nstcr = 0
    evalsp = Vector{Vector{Float64}}(undef,Nspecies)
    occsp  = Vector{Vector{Float64}}(undef,Nspecies)

    # XXX Will be setup properly later
    rsp = Vector{Vector{Float64}}(undef,Nspecies)
    rhosp = Vector{Vector{Float64}}(undef,Nspecies)
    vrsp = Vector{Vector{Float64}}(undef,Nspecies)

    xctsp  = (3,0,0) # not used, always defaulting to this?

    return AtomicSpeciesVars(
        sppath, spfname, spname, spsymb, spzn,
        ptnucl, rnucl, volnucl, nrnucl, nrcnucl,
        vcln, spze, spmass, rminsp, rmaxsp, nrsp,
        maxstsp, nstsp, ecvcut, esccut,
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