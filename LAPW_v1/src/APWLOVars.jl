mutable struct APWLOVars
    # energy step used for numerical calculation of energy derivatives
    deapwlo::Float64
    # maximum allowable APW order
    maxapword::Int64 # parameter=4
    # APW order
    apword::Array{Int64,2} #(0:maxlapw,maxspecies)
    # maximum of apword over all angular momenta and species
    apwordmax::Int64
    # total number of APW coefficients (l, m and order) for each species
    lmoapw::Vector{Int64} # (maxspecies)
    # polynomial order used for APW radial derivatives
    npapw::Int64
    # APW initial linearisation energies
    apwe0::Array{Float64,3} #(maxapword,0:maxlapw,maxspecies)
    # APW linearisation energies
    apwe::Array{Float64,3} # (:,:,:)
    # APW derivative order
    apwdm::Array{Int64,3} # (maxapword,0:maxlapw,maxspecies)
    # apwve is .true. if the linearisation energies are allowed to vary
    apwve::Array{Bool,3} # (maxapword,0:maxlapw,maxspecies)
    # APW radial functions
    apwfr::Array{Float64,5}
    # derivate of radial functions at the muffin-tin surface
    apwdfr::Array{Float64,3}
    # maximum number of local-orbitals
    maxlorb::Int64 # parameter=200
    # maximum allowable local-orbital order
    maxlorbord::Int64 # parameter=5
    # number of local-orbitals
    nlorb::Vector{Int64} #(maxspecies)
    # maximum nlorb over all species
    nlomax::Int64
    # total number of local-orbitals
    nlotot::Int64
    # local-orbital order
    lorbord::Array{Int64,2} #(maxlorb,maxspecies)
    # maximum lorbord over all species
    lorbordmax::Int64
    # polynomial order used for local-orbital radial derivatives
    nplorb::Int64
    # local-orbital angular momentum
    lorbl::Array{Float64,2} #(maxlorb,maxspecies)
    # maximum lorbl over all species
    lolmax::Int64
    # (lolmax+1)^2
    lolmmax::Int64
    # local-orbital initial energies
    lorbe0::Array{Float64,3} # (maxlorbord,maxlorb,maxspecies)
    # local-orbital energies
    lorbe::Array{Float64,3} #(:,:,:)
    # local-orbital derivative order
    lorbdm::Array{Int64,3} #(maxlorbord,maxlorb,maxspecies)
    # lorbve is .true. if the linearisation energies are allowed to vary
    lorbve::Array{Bool,3} # (maxlorbord,maxlorb,maxspecies)
    # local-orbital radial functions
    lofr::Array{Float64,4} # (:,:,:,:)
    # band energy search tolerance
    epsband::Float64
    # maximum allowed change in energy during band energy search; enforced only if
    # default energy is less than zero
    demaxbnd::Float64
    # minimum default linearisation energy over all APWs and local-orbitals
    e0min::Float64
    # if autolinengy is .true. then the fixed linearisation energies are set to the
    # Fermi energy minus dlefe
    autolinengy::Bool
    # difference between linearisation and Fermi energies when autolinengy is .true.
    dlefe::Float64
    # lorbcnd is .true. if conduction state local-orbitals should be added
    lorbcnd::Bool
    # conduction state local-orbital order
    lorbordc::Bool
    # excess order of the APW and local-orbital functions
    nxoapwlo::Int64
    # excess local orbitals
    nxlo::Int64
end

function debug_apwlo(
    atoms, specs_info::Vector{SpeciesInfo}, mt_vars::MuffinTins
)
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    # The actual lmax used in APW
    lmaxapw = mt_vars.lmaxapw

    # copy from specs_info
    #
    #
    # This array depend on angular momentum (?)
    apword = Vector{OffsetVector{Int64,Vector{Int64}}}(undef, Nspecies)
    for isp in 1:Nspecies
        apword[isp] = OffsetArray(zeros(Int64,lmaxapw+1), 0:lmaxapw) # allocate
        apword[isp][:] = specs_info[isp].apword[0:lmaxapw] # copy
    end
    #
    # These arrays depend on order and angular momentum index (?)
    #
    apwe0 = Vector{OffsetMatrix{Float64,Matrix{Float64}}}(undef, Nspecies)
    for isp in 1:Nspecies
        maxapword = specs_info[isp].maxapword
        apwe0[isp] = OffsetArray(zeros(Float64, maxapword, lmaxapw+1), 1:maxapword, 0:lmaxapw)
        apwe0[isp][:] = specs_info[isp].apwe0[1:maxapword,0:lmaxapw]
    end
    #
    apwdm = Vector{OffsetMatrix{Int64,Matrix{Int64}}}(undef, Nspecies)
    for isp in 1:Nspecies
        maxapword = specs_info[isp].maxapword
        apwdm[isp] = OffsetArray(zeros(Int64, maxapword, lmaxapw+1), 1:maxapword, 0:lmaxapw)
        apwdm[isp][:] = specs_info[isp].apwdm[1:maxapword,0:lmaxapw]
    end
    apwve = Vector{OffsetMatrix{Bool,Matrix{Bool}}}(undef, Nspecies)
    for isp in 1:Nspecies
        maxapword = specs_info[isp].maxapword
        apwve[isp] = OffsetArray(zeros(Bool, maxapword, lmaxapw+1), 1:maxapword, 0:lmaxapw)
        apwve[isp][:] = specs_info[isp].apwve[1:maxapword,0:lmaxapw]
    end
    #
    nlorb = zeros(Int64, Nspecies)
    for isp in 1:Nspecies
        nlorb[isp] = specs_info[isp].nlorb
    end
    #
    lorbl = Vector{Vector{Int64}}(undef, Nspecies)
    lorbord = Vector{Vector{Int64}}(undef, Nspecies)
    for isp in 1:Nspecies
        # Preallocate arrays?
        lorbl[isp] = deepcopy(specs_info[isp].lorbl)
        lorbord[isp] = deepcopy(specs_info[isp].lorbord)
    end
    #
    lorbe0 = Vector{Vector{Vector{Float64}}}(undef, Nspecies)
    lorbdm = Vector{Vector{Vector{Int64}}}(undef, Nspecies)
    lorbve = Vector{Vector{Vector{Bool}}}(undef, Nspecies)
    for isp in 1:Nspecies
        # XXX explicit copy?
        lorbe0[isp] = deepcopy(specs_info[isp].lorbe0)
        lorbdm[isp] = deepcopy(specs_info[isp].lorbdm)
        lorbve[isp] = deepcopy(specs_info[isp].lorbve)
    end

    lmoapw = zeros(Int, Nspecies)
    apwordmax = 0
    lorbordmax = 0
    nlomax = 0
    lolmax = 0
    for isp in 1:Nspecies
        lmoapw[isp] = 0
        for l1 in 0:lmaxapw
            # find the maximum APW order
            apwordmax = max(apwordmax, apword[isp][l1])
            # find total number of APW coefficients (l, m and order)
            lmoapw[isp] += (2*l1 + 1)*apword[isp][l1]
        end
        # find the maximum number of local-orbitals
        nlomax = max(nlomax, nlorb[isp])
        # find the maximum local-orbital order and angular momentum
        for ilo in 1:nlorb[isp]
            lolmax = max( lolmax, lorbl[isp][ilo] )
            lorbordmax = max( lorbordmax, lorbord[isp][ilo] )
        end
    end
    lolmmax = (lolmax + 1)^2

    # polynomial order used for APW and local-orbital radial derivatives
    npapw = max(apwordmax+1, 4)
    nplorb = max(lorbordmax+1, 4)

    # set the APW and local-orbital linearisation energies to the default
    apwe = Vector{OffsetMatrix{Float64,Matrix{Float64}}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        maxapword = specs_info[isp].maxapword
        apwe[ia] = OffsetArray(zeros(Float64, maxapword, lmaxapw+1), 1:maxapword, 0:lmaxapw)
        for l1 in 0:lmaxapw, io in 1:apword[isp][l1]
            apwe[ia][io,l1] = apwe0[isp][io,l1]
        end
    end

    lorbe = Vector{Vector{Vector{Float64}}}(undef, Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        lorbe[ia] = similar.(lorbe0[isp]) # should be safe for Vector
        for ilo in 1:nlorb[isp], io in 1:lorbord[isp][ilo]
            lorbe[ia][ilo][io] = lorbe0[isp][ilo][io]
        end
    end

    # wavefunctions ....


    @infiltrate

    return
end



# effective size for 0:maxlapw -> 0:lmaxapw
# Using indices starting from 1
function APWLOVars(Nspecies::Int64, maxlapw::Int64)

    deapwlo = 0
    maxapword = 4
    apword = zeros(Int64,maxlapw+1,Nspecies)
    apwordmax = 0
    lmoapw = zeros(Int64,Nspecies)
    npapw = 0
    
    apwe0 = zeros(Float64, maxapword, maxlapw+1, Nspecies)
    apwe = zeros(Float64, 1, 1, 1)
    apwdm = zeros(Int64,maxapword, maxlapw+1, Nspecies)
    
    apwve = zeros(Bool, maxapword, maxlapw+1, Nspecies)
    apwfr = zeros(Float64,1,1,1,1,1)
    apwdfr = zeros(Float64,1,1,1)

    maxlorb = 200
    maxlorbord = 5
    nlorb = zeros(Int64,Nspecies)
    nlomax = 0
    nlotot = 0
    lorbord = zeros(Int64,maxlorb,Nspecies)
    lorbordmax = 0
    nplorb = 0

    lorbl = zeros(Float64,maxlorb,Nspecies)
    lolmax = 0
    lolmmax = 0
    lorbe0 = zeros(Float64,maxlorbord,maxlorb,Nspecies)
    lorbe = zeros(Float64,1,1,1)
    lorbdm = zeros(Int64,maxlorbord,maxlorb,Nspecies)
    lorbve = zeros(Bool,maxlorbord,maxlorb,Nspecies)
    lofr = zeros(Float64,1,1,1,1)
    epsband = 0.0
    demaxbnd = 0.0
    e0min = 0.0
    autolinengy = false 
    dlefe = 0.0
    lorbcnd = false
    lorbordc = 0
    nxoapwlo = 0
    nxlo = 0

    return APWLOVars(
        deapwlo, maxapword, apword, apwordmax,
        lmoapw, npapw, apwe0, apwe, apwdm, apwve, apwfr, apwdfr,
        maxlorb, maxlorbord, nlorb, nlomax, nlotot, lorbord,
        lorbordmax, nplorb, lorbl, lolmax,
        lolmmax, lorbe0, lorbe, lorbdm, lorbve,
        lofr, epsband, demaxbnd, e0min, autolinengy,
        dlefe, lorbcnd, lorbordc, nxoapwlo, nxlo
    )

end

#apwlo_vars = APWLOVars(2,50)

# indexing for array allocated with lmaxo
# allocate array((-lmaxo-1:lmaxo+2))
# from lmaxo idx to the usual 1-based indexing
function lmaxo2idx(idx, lmaxo)
    return idx + lmaxo + 2
end