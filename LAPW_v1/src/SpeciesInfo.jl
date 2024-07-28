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
    apwe0::OffsetMatrix{Float64,Matrix{Float64}}
    apwdm::OffsetMatrix{Int64,Matrix{Int64}}
    apwve::OffsetMatrix{Bool,Matrix{Bool}}
    nlx::Int64
    nlorb::Int64
    lorbl::Vector{Int64}
    lorbord::Vector{Int64}
    lorbe0::Vector{Vector{Float64}}
    lorbdm::Vector{Vector{Int64}}
    lorbve::Vector{Vector{Bool}}
end



function SpeciesInfo(
    filename::String;
    verbose=true,
    maxlapw=50,
    maxapword=4,
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
    apword[1:maxlapw] .= apword[0]
    # XXX Generally, apword will depend on l ? Here, the order is the same for all l.
    #
    e0min = 0.0
    apwe0 = OffsetArray( zeros(Float64, maxapword, maxlapw+1), 1:maxapword, 0:maxlapw )
    apwdm = OffsetArray( zeros(Int64, maxapword, maxlapw+1), 1:maxapword, 0:maxlapw )
    apwve = OffsetArray( zeros(Bool, maxapword, maxlapw+1), 1:maxapword, 0:maxlapw )
    # apwe0, apwdm, apwve
    #XXX In all species files considered this is only one line
    for iord in 1:apword[0]
        line = readline(f)
        ll = split(line)
        apwe0[iord,0] = parse(Float64, ll[1])
        apwdm[iord,0] = parse(Int64, ll[2])
        if ll[3] == "T"
            apwve[iord,0] = true
        elseif ll[3] == "F"
            apwve[iord,0] = false
        else
            error("Unable to parse apwve")
        end
        # 
        # set the APW linearisation energies, derivative orders and variability for l > 0
        apwe0[iord,1:maxlapw] .= apwe0[iord,0]
        apwdm[iord,1:maxlapw] .= apwdm[iord,0]
        apwve[iord,1:maxlapw] .= apwve[iord,0]
        #
        e0min = min( e0min, apwe0[iord,0] )
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
        end
    end

    close(f)

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