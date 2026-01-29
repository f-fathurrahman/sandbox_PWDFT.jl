mutable struct SpeciesInfo
    #
    filename::String # usually contains .in prefix
    spsymb::String # the actual symbol, e.g. spsymb
    spname::String # full name, e.g. silicon, platinum, etc
    spzn::Float64
    spmass::Float64
    rminsp::Float64
    rmaxsp::Float64
    #
    rmt::Float64
    nrmt::Int64
    #
    nstsp::Int64
    nsp::Vector{Int64}
    lsp::Vector{Int64}
    ksp::Vector{Int64}
    spcore::Vector{Bool}
    occsp::Vector{Float64}
    evalsp::Vector{Float64}
    #
    maxlapw::Int64
    maxapword::Int64
    maxlorb::Int64
    maxlorbord::Int64
    #
    e0min::Float64
    apword::OffsetVector{Int64,Vector{Int64}}
    apwe0::OffsetVector{Vector{Float64}, Vector{Vector{Float64}}}
    apwdm::OffsetVector{Vector{Int64}, Vector{Vector{Int64}}}
    apwve::OffsetVector{Vector{Bool}, Vector{Vector{Bool}}}
    nlx::Int64
    nlorb::Int64
    lorbl::Vector{Int64}
    lorbord::Vector{Int64}
    lorbe0::Vector{Vector{Float64}}
    lorbdm::Vector{Vector{Int64}}
    lorbve::Vector{Vector{Bool}}
end



function SpeciesInfo(
    filename;
    verbose=true,
    maxlapw=50,
    maxlorb=200,
    maxlorbord=5
)

    f = open(filename, "r")
    
    # species symbol
    line = readline(f)
    spsymb = replace(split(line)[1], "'" => "")
    verbose && println("spsymb = ", spsymb)
    
    # species name
    line = readline(f)
    spname = replace(split(line)[1], "'" => "")
    verbose && println("spname = ", spname)
    
    # atomic number
    line = readline(f)
    spzn = parse(Float64, split(line)[1])
    verbose && println("spzn   = ", spzn)

    # mass
    line = readline(f)
    spmass = parse(Float64, split(line)[1])
    verbose && println("spmass = ", spmass)

    # Radial mesh
    line = readline(f)
    ll = split(line)
    
    rminsp = parse(Float64, ll[1])
    
    rmt = parse(Float64, ll[2])  # muffin tin
    
    rmaxsp = parse(Float64, ll[3])
    
    nrmt = parse(Int64, ll[4])
    
    verbose && println("rminsp = ", rminsp)
    verbose && println("rmt    = ", rmt)
    verbose && println("rmaxsp = ", rmaxsp)
    verbose && println("nrmt   = ", nrmt)

    # Atomic states
    line = readline(f)
    nstsp = parse(Int64, split(line)[1])    
    nsp = zeros(Int64,nstsp)
    lsp = zeros(Int64,nstsp)
    ksp = zeros(Int64,nstsp)
    spcore = zeros(Bool,nstsp)
    occsp = zeros(Float64,nstsp)
    evalsp = zeros(Float64,nstsp)
    #
    for ist in 1:nstsp
        #
        line = readline(f)
        ll = split(line)
        #
        nsp[ist] = parse(Int64, ll[1])
        lsp[ist] = parse(Int64, ll[2])
        ksp[ist] = parse(Int64, ll[3]) # whats this?
        #
        occsp[ist] = parse(Float64, ll[4])
        if ll[5] == "T"
            spcore[ist] = true
        elseif ll[5] == "F"
            spcore[ist] = false
        else
            error("Unable to parse spcore")
        end
    end

    # XXX: maxlapw=50 may be too large
    # maxlapw is a hardcoded parameter in m_muffin_tins, used to hardcoded some array size
    # The actual value is lmaxapw which is set to 8 in default_muffin_tins

    # apword
    line = readline(f)
    apword = OffsetArray( zeros(Int64, maxlapw+1), 0:maxlapw )
    apword[0] = parse(Int64, split(line)[1])
    # set the APW orders for > 0
    apword[1:maxlapw] .= apword[0] # XXX apword are the same for all l ?
    #
    maxapword = maximum(apword)
    MAX_APW_ORD = 4 # Hardcoded, probably not needed
    @assert maxapword <= MAX_APW_ORD
    #
    # XXX Generally, apword will depend on l ? Here, the order is the same for all l.
    #
    e0min = 0.0
    apwe0 = OffsetArray( Vector{Vector{Float64}}(undef, maxlapw+1), 0:maxlapw )
    apwdm = OffsetArray( Vector{Vector{Int64}}(undef, maxlapw+1), 0:maxlapw )
    apwve = OffsetArray( Vector{Vector{Bool}}(undef, maxlapw+1), 0:maxlapw )
    # apwe0, apwdm, apwve
    for l in 0:maxlapw
        apwe0[l] = zeros(Float64, apword[l])
        apwdm[l] = zeros(Int64, apword[l])
        apwve[l] = zeros(Bool, apword[l])
    end
    #XXX In all species files considered this is only one line
    for iord in 1:apword[0] # only apword[0] ?
        line = readline(f)
        ll = split(line)
        apwe0[0][iord] = parse(Float64, ll[1])
        apwdm[0][iord] = parse(Int64, ll[2])
        if ll[3] == "T"
            apwve[0][iord] = true
        elseif ll[3] == "F"
            apwve[0][iord] = false
        else
            error("Unable to parse apwve")
        end
        # 
        # set the APW linearization energies, derivative orders and variability for l > 0
        for l in 1:maxlapw
            apwdm[l][iord] = apwdm[0][iord]
            apwe0[l][iord] = apwe0[0][iord]
            apwve[l][iord] = apwve[0][iord]
        end
        e0min = min( e0min, apwe0[0][iord] )
    end

    # nlx (skipped?)
    line = readline(f)
    nlx = parse(Int64, split(line)[1])
    if nlx > 0
        error("Unsupported value of nlx")
    end

    # nlorb
    line = readline(f)
    nlorb = parse(Int64, split(line)[1])
    @assert nlorb <= maxlorb
    #
    lorbl = zeros(Int64, nlorb)
    lorbord = zeros(Int64, nlorb)
    lorbe0 = Vector{Vector{Float64}}(undef, nlorb)
    lorbdm = Vector{Vector{Int64}}(undef, nlorb)
    lorbve = Vector{Vector{Bool}}(undef, nlorb)
    for iorb in 1:nlorb
        # lorbl, lorbord
        line = readline(f); ll = split(line)
        #
        lorbl[iorb] = parse(Int64, ll[1])
        #
        lorbord[iorb] = parse(Int64, ll[2])
        @assert lorbord[iorb] <= maxlorbord
        #
        lorbe0[iorb] = zeros(Float64, lorbord[iorb])
        lorbdm[iorb] = zeros(Int64, lorbord[iorb])
        lorbve[iorb] = zeros(Bool, lorbord[iorb])
        #
        for i in 1:lorbord[iorb]
            # lorbe0, lorbdm, lorbve
            line = readline(f); ll = split(line)
            #
            lorbe0[iorb][i] = parse(Float64, ll[1])
            #
            lorbdm[iorb][i] = parse(Float64, ll[2])
            #
            if ll[3] == "T"
                lorbve[iorb][i] = true
            elseif ll[3] == "F"
                lorbve[iorb][i] = false
            else
                error("Unable to parse lorbve")
            end
            #
            e0min = min(e0min, lorbe0[iorb][i])
        end
    end

    close(f)

    # subtract 2 Hartree from the minimum energy
    e0min = e0min - 2.0

    return SpeciesInfo(
        filename, spsymb, spname, spzn, spmass,
        rminsp, rmaxsp,
        rmt, nrmt,
        nstsp, nsp, lsp, ksp, spcore, occsp, evalsp,
        maxlapw, maxapword, maxlorb, maxlorbord,
        e0min, apword, apwe0, apwdm, apwve,
        nlx, nlorb,
        lorbl, lorbord, lorbe0, lorbdm, lorbve
    )

end