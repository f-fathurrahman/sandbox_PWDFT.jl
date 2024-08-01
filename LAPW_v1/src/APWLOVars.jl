mutable struct APWLOVars
    #
    # energy step used for numerical calculation of energy derivatives
    deapwlo::Float64
    #
    # APW order
    apword::Vector{OffsetVector{Int64, Vector{Int64}}}
    #
    # maximum of apword over all angular momenta and species
    apwordmax::Int64
    #
    # total number of APW coefficients (l, m and order) for each species
    lmoapw::Vector{Int64}
    #
    # polynomial order used for APW radial derivatives
    npapw::Int64
    #
    # APW initial linearisation energies
    apwe0::Vector{OffsetVector{Vector{Float64}, Vector{Vector{Float64}}}}
    #
    # APW linearisation energies
    apwe::Vector{OffsetVector{Vector{Float64}, Vector{Vector{Float64}}}}
    #
    # APW derivative order
    apwdm::Vector{OffsetVector{Vector{Int64}, Vector{Vector{Int64}}}}
    #
    # apwve is .true. if the linearisation energies are allowed to vary
    apwve::Vector{OffsetVector{Vector{Bool}, Vector{Vector{Bool}}}}
    #
    # APW radial functions
    apwfr::Vector{OffsetMatrix{Matrix{Float64}, Matrix{Matrix{Float64}}}}
    #
    # derivate of radial functions at the muffin-tin surface
    apwdfr::Vector{OffsetMatrix{Float64, Matrix{Float64}}}
    #
    # maximum number of local-orbitals
    maxlorb::Int64 # parameter=200
    #
    # maximum allowable local-orbital order
    maxlorbord::Int64 # parameter=5
    #
    # number of local-orbitals
    nlorb::Vector{Int64} #(maxspecies)
    #
    # maximum nlorb over all species
    nlomax::Int64
    #
    # total number of local-orbitals
    nlotot::Int64
    # local-orbital order
    lorbord::Vector{Vector{Int64}}
    #
    # maximum lorbord over all species
    lorbordmax::Int64
    #
    # polynomial order used for local-orbital radial derivatives
    nplorb::Int64
    #
    # local-orbital angular momentum
    lorbl::Vector{Vector{Int64}}
    #
    # maximum lorbl over all species
    lolmax::Int64
    #
    # (lolmax+1)^2
    lolmmax::Int64
    # local-orbital initial energies
    lorbe0::Vector{Vector{Vector{Float64}}}
    #
    # local-orbital energies
    lorbe::Vector{Vector{Vector{Float64}}}
    #
    # local-orbital derivative order
    lorbdm::Vector{Vector{Vector{Int64}}}
    #
    # lorbve is .true. if the linearization energies are allowed to vary
    lorbve::Vector{Vector{Vector{Bool}}}
    #
    # local-orbital radial functions
    lofr::Vector{Vector{Matrix{Float64}}}
    #
    # band energy search tolerance
    epsband::Float64
    #
    # maximum allowed change in energy during band energy search; enforced only if
    # default energy is less than zero
    demaxbnd::Float64
    #
    # minimum default linearisation energy over all APWs and local-orbitals
    e0min::Float64
    #
    # if autolinengy is .true. then the fixed linearisation energies are set to the
    # Fermi energy minus dlefe
    autolinengy::Bool
    #
    # difference between linearisation and Fermi energies when autolinengy is .true.
    dlefe::Float64
    #
    # lorbcnd is .true. if conduction state local-orbitals should be added
    lorbcnd::Bool
    #
    # conduction state local-orbital order
    lorbordc::Int64
    #
    # excess order of the APW and local-orbital functions
    nxoapwlo::Int64
    #
    # excess local orbitals
    nxlo::Int64
end

function APWLOVars(
    atoms, specs_info::Vector{SpeciesInfo}, mt_vars::MuffinTins;
    deapwlo = 0.05,
    maxlorb = 200,
    maxlorbord = 5,
    epsband = 1e-12,
    demaxbnd = 2.5,
    autolinengy = false,
    dlefe = -0.1,
    lorbcnd = false,
    lorbordc = 3,
    nxoapwlo = 0,
    nxlo = 0
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
    APWE0_ELTYPE = typeof(specs_info[1].apwe0)
    apwe0 = Vector{APWE0_ELTYPE}(undef, Nspecies)
    for isp in 1:Nspecies
        apwe0[isp] = OffsetArray( Vector{Vector{Float64}}(undef, lmaxapw+1), 0:lmaxapw )
        for l in 0:lmaxapw
            apwe0[isp][l] = zeros(Float64, apword[isp][l])
            for io in 1:apword[isp][l]
                apwe0[isp][l][io] = specs_info[isp].apwe0[l][io]
            end
        end
    end
    #
    APWDM_ELTYPE = typeof(specs_info[1].apwdm)
    apwdm = Vector{APWDM_ELTYPE}(undef, Nspecies)
    for isp in 1:Nspecies
        apwdm[isp] = OffsetArray( Vector{Vector{Int64}}(undef, lmaxapw+1), 0:lmaxapw )
        for l in 0:lmaxapw
            apwdm[isp][l] = zeros(Float64, apword[isp][l])
            for io in 1:apword[isp][l]
                apwdm[isp][l][io] = specs_info[isp].apwdm[l][io]
            end
        end
    end
    APWVE_ELTYPE = typeof(specs_info[1].apwve)
    apwve = Vector{APWVE_ELTYPE}(undef, Nspecies)
    for isp in 1:Nspecies
        apwve[isp] = OffsetArray( Vector{Vector{Bool}}(undef, lmaxapw+1), 0:lmaxapw )
        for l in 0:lmaxapw
            apwve[isp][l] = zeros(Float64, apword[isp][l])
            for io in 1:apword[isp][l]
                apwve[isp][l][io] = specs_info[isp].apwve[l][io]
            end
        end
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
    APWE_ELTYPE = eltype(apwe0)
    apwe = Vector{APWE_ELTYPE}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        apwe[isp] = OffsetArray( Vector{Vector{Float64}}(undef, lmaxapw+1), 0:lmaxapw )
        for l in 0:lmaxapw
            apwe[ia][l] = zeros(Float64, apword[isp][l])
            for io in 1:apword[isp][l]
                apwe[ia][l][io] = apwe0[isp][l][io]
            end
        end
    end
    #
    lorbe = Vector{Vector{Vector{Float64}}}(undef, Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        lorbe[ia] = similar.(lorbe0[isp]) # should be safe for Vector
        for ilo in 1:nlorb[isp], io in 1:lorbord[isp][ilo]
            lorbe[ia][ilo][io] = lorbe0[isp][ilo][io]
        end
        # ilo: local orbital index
        # io: order index
    end

    #
    # allocate radial function arrays
    #
    #XXX Ugh... This type is complicated
    #XXX Probably better use simple multidimensional array
    APW_RADIAL_TYPE = OffsetMatrix{Matrix{Float64}, Matrix{Matrix{Float64}}}
    # An OffsetMatrix whose element is a matrix
    #
    apwfr = Vector{APW_RADIAL_TYPE}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = mt_vars.nrmt[isp]
        maxapword = specs_info[isp].maxapword
        apwfr[ia] = OffsetArray(
            Array{Matrix{Float64}}(undef, maxapword, lmaxapw+1), 1:maxapword, 0:lmaxapw
        )
        for l1 in 0:lmaxapw, io in 1:apword[isp][l1]
            apwfr[ia][io,l1] = zeros(Float64, nr, 2)
        end
    end
    # apwfr[atom index][order index, ang mom index][radial mt index,major-minor?]
    # XXX: order index is singleton anyway
    apwdfr = Vector{OffsetMatrix{Float64, Matrix{Float64}}}(undef,Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        maxapword = specs_info[isp].maxapword
        apwdfr[ia] = OffsetArray(zeros(Float64, maxapword, lmaxapw+1), 1:maxapword, 0:lmaxapw)
    end

    lofr = Vector{Vector{Matrix{Float64}}}(undef, Natoms)
    for ia in 1:Natoms
        isp = atm2species[ia]
        nr = mt_vars.nrmt[isp]
        lofr[ia] = Vector{Matrix{Float64}}(undef,nlorb[isp])
        for ilo in 1:nlorb[isp]
            lofr[ia][ilo] = zeros(Float64,nr,2)
        end
    end

    nlotot = 0
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ilo in 1:nlorb[isp]
            l = lorbl[isp][ilo]
            for m in -l:l
                nlotot += 1
                #lm = mt_vars.idxlm[l,m]
                #idxlo[lm,ilo,ias] = nlotot
            end
        end
    end
    @info "nlotot = $(nlotot)"

    e0min = 0.0
    for isp in 1:Nspecies
        if specs_info[isp].e0min < e0min
            e0min = specs_info[isp].e0min
        end
    end
    @info "e0min = $(e0min)"

    #@infiltrate

    return APWLOVars(
        deapwlo, apword, apwordmax,
        lmoapw, npapw, apwe0, apwe, apwdm, apwve, apwfr, apwdfr,
        maxlorb, maxlorbord, nlorb, nlomax, nlotot, lorbord,
        lorbordmax, nplorb, lorbl, lolmax,
        lolmmax, lorbe0, lorbe, lorbdm, lorbve,
        lofr, epsband, demaxbnd, e0min, autolinengy,
        dlefe, lorbcnd, lorbordc, nxoapwlo, nxlo
    )
end

