# NOTES:
# Several variables are array of size Nspecies (specific for each species)
# Other variables are common for both species

mutable struct MuffinTins
    # scale factor for number of muffin-tin points
    nrmtscf::Float64
    # number of muffin-tin radial points for each species
    nrmt::Vector{Int64}
    # optional default muffin-tin radius for all atoms
    rmtall::Float64
    # minimum allowed distance between muffin-tin surfaces
    rmtdelta::Float64
    # muffin-tin radii
    rmt::Vector{Float64}
    # total muffin-tin volume
    omegamt::Float64
    # radial step length for coarse mesh
    lradstp::Int64
    # number of coarse radial mesh points
    nrcmt::Vector{Int64}
    # coarse muffin-tin radial mesh
    rcmt::Vector{Vector{Float64}}
    # r^l on fine radial mesh
    rlmt::Vector{OffsetArray{Float64,2,Array{Float64,2}}}
    # r^l on coarse radial mesh
    rlcmt::Vector{OffsetArray{Float64,2,Array{Float64,2}}}
    # weights for spline integration on fine radial mesh
    wrmt::Vector{Vector{Float64}}
    # weights for spline partial integration on fine radial mesh
    wprmt::Vector{Matrix{Float64}}
    # weights for spline integration on coarse radial mesh
    wrcmt::Vector{Vector{Float64}}
    # weights for spline partial integration on coarse radial mesh
    wprcmt::Vector{Matrix{Float64}}
    # maximum allowable angular momentum for augmented plane waves
    maxlapw::Int64 # parameter =50
    # maximum angular momentum for augmented plane waves
    lmaxapw::Int64
    # (lmaxapw+1)^2
    lmmaxapw::Int64
    # maximum angular momentum on the outer part of the muffin-tin
    lmaxo::Int64
    # (lmaxo+1)^2
    lmmaxo::Int64
    # maximum angular momentum on the inner part of the muffin-tin
    lmaxi::Int64
    # (lmaxi+1)^2
    lmmaxi::Int64
    # fraction of muffin-tin radius which constitutes the inner part
    fracinr::Float64
    # number of fine/coarse radial points on the inner part of the muffin-tin
    nrmti::Vector{Int64}
    nrcmti::Vector{Int64} 
    # index to (l,m) pairs
    idxlm::OffsetMatrix{Int64, Matrix{Int64}}
    # inverse index to (l,m) pairs
    idxil::Vector{Int64}
    idxim::Vector{Int64}
    # number of fine/coarse points in packed muffin-tins
    npmti::Vector{Int64}
    npmt::Vector{Int64}
    npcmti::Vector{Int64}
    npcmt::Vector{Int64}
    #
    SHT::SphericalHarmonicTransform
    #
    gntyry::Array{ComplexF64,3}
end


function MuffinTins(
    specs_info::Vector{SpeciesInfo}, atsp_vars;
    lmaxi=1,
    rmtdelta=0.05,
    lradstp=4,
    maxlapw=50,
    lmaxapw=8,
    lmaxo=6,
    fracinr=0.01
)
    @assert lmaxo <= lmaxapw

    Nspecies = length(specs_info)

    nrmt = zeros(Int64, Nspecies)
    nrcmt = zeros(Int64, Nspecies)
    rmt = zeros(Float64, Nspecies)
    for isp in 1:Nspecies
        nrmt[isp] = specs_info[isp].nrmt
        @assert mod(nrmt[isp]-1, lradstp) == 0 # should be commensurate
        #
        nrcmt[isp] = (nrmt[isp] - 1)/lradstp + 1
        #
        rmt[isp] = specs_info[isp].rmt
    end

    lmmaxapw = (lmaxapw+1)^2
    lmmaxo = (lmaxo+1)^2
    lmaxi  = min(lmaxi,lmaxo)
    lmmaxi = (lmaxi+1)^2
    # index to (l,m) pairs
    idxlm = OffsetArray( zeros(Int64,lmaxapw+1,2*lmaxapw+1), 0:lmaxapw, -lmaxapw:lmaxapw)
    idxil = zeros(Int64, lmmaxapw)
    idxim = zeros(Int64, lmmaxapw)
    lm = 0
    for l in 0:lmaxapw, m in -l:l
        lm = lm + 1
        idxlm[l,m] = lm
        idxil[lm] = l
        idxim[lm] = m
    end


    rlmt = Vector{OffsetMatrix{Float64, Matrix{Float64}}}(undef,Nspecies)
    wrmt = Vector{Vector{Float64}}(undef,Nspecies)
    wprmt = Vector{Matrix{Float64}}(undef,Nspecies)
    for isp in 1:Nspecies
        rlmt[isp] = OffsetArray(
            zeros(Float64, nrmt[isp], 2*(lmaxo+2)),
            1:nrmt[isp], -lmaxo-1:lmaxo+2
        )
        wrmt[isp] = zeros(Float64, nrmt[isp])
        wprmt[isp] = zeros(Float64, 4, nrmt[isp])
    end
    # Now fill in the values
    rsp = atsp_vars.rsp
    #
    for isp in 1:Nspecies
        #
        # calculate r^l on the fine radial mesh
        nr = nrmt[isp]
        @assert nr <= atsp_vars.nrsp[isp] # XXX Need this?
        #
        for ir in 1:nr
            rlmt[isp][ir,-1] = 1.0/rsp[isp][ir]
            rlmt[isp][ir,0] = 1.0
            rlmt[isp][ir,1] = rsp[isp][ir]
        end
        #
        for l in range(-2,stop=-lmaxo-1,step=-1)
            for ir in 1:nr
                rlmt[isp][ir,l] = rlmt[isp][ir,l+1]/rsp[isp][ir]
            end
        end
        #
        for l in 2:lmaxo+2            
            for ir in 1:nr
                rlmt[isp][ir,l] = rlmt[isp][ir,l-1] * rsp[isp][ir]
            end
        end
        # determine the weights for spline integration on the fine radial mesh
        wsplint!(nr, rsp[isp], wrmt[isp])
        # multiply by r^2
        for ir in 1:nr
            wrmt[isp][ir] = wrmt[isp][ir]*rlmt[isp][ir,2]
        end
        # determine the weights for partial integration on fine radial mesh
        wsplintp!(nr, rsp[isp], wprmt[isp])
    end


    # determine the fraction of the muffin-tin radius which defines the inner part
    if fracinr < 0.0
        # be explicit about conversion to Float64 (not really needed actually for Julia)
        fracinr = sqrt( Float64(lmmaxi) / Float64(lmmaxo) )
    end


    nrmtscf = 0
    rmtall = 0
    omegamt = 0.0

    # Coarse mesh
    rcmt = Vector{Vector{Float64}}(undef,Nspecies)
    rlcmt = Vector{OffsetArray{Float64,2,Array{Float64,2}}}(undef,Nspecies)
    wrcmt = Vector{Vector{Float64}}(undef,Nspecies)
    wprcmt = Vector{Matrix{Float64}}(undef,Nspecies)
    # set up the coarse radial meshes and find the inner part of the muffin-tin
    # where rho is calculated with lmaxi
    for isp in 1:Nspecies
        rcmt[isp] = zeros(Float64, nrcmt[isp])
        rlcmt[isp] = OffsetArray(
            zeros(Float64, nrcmt[isp], 2*(lmaxo+2)),
            1:nrcmt[isp], -lmaxo-1:lmaxo+2
        )
        wrcmt[isp] = zeros(Float64, nrcmt[isp])
        wprcmt[isp] = zeros(Float64, 4, nrcmt[isp])
    end
    #
    # Inner radial grid
    nrmti = zeros(Int64, Nspecies)
    nrcmti = zeros(Int64, Nspecies)
    #
    for isp in 1:Nspecies
        t1 = fracinr*rmt[isp]
        nrmti[isp] = 1
        nrcmti[isp] = 1 # need this?
        irc = 0
        for ir in range(1,stop=nrmt[isp],step=lradstp)
            irc = irc + 1
            rcmt[isp][irc] = rsp[isp][ir]
            if rsp[isp][ir] < t1
                nrmti[isp] = ir
                nrcmti[isp] = irc
            end
        end
        # store r^l on the coarse radial mesh
        for l in range(-lmaxo-1, stop=lmaxo+2)
            irc = 0
            for ir in range(1,stop=nrmt[isp],step=lradstp)
                irc = irc + 1
                rlcmt[isp][irc,l] = rlmt[isp][ir,l]
            end
        end
        # determine the weights for spline integration on the coarse radial mesh
        nrc = nrcmt[isp]
        wsplint!(nrc, rcmt[isp], wrcmt[isp])
        # multiply by r^2
        for ir in 1:nrc
            wrcmt[isp][ir] = wrcmt[isp][ir] * rlcmt[isp][ir,2]
        end
        # determine the weights for partial integration on coarse radial mesh
        wsplintp!(nrc, rcmt[isp], wprcmt[isp])
    end


    # packed MT
    npmti = zeros(Int64, Nspecies)
    npmt = zeros(Int64, Nspecies)
    npcmti = zeros(Int64, Nspecies)
    npcmt = zeros(Int64, Nspecies)
    for isp in 1:Nspecies
        #
        npmti[isp] = lmmaxi*nrmti[isp]
        npmt[isp] = npmti[isp] + lmmaxo*(nrmt[isp] - nrmti[isp])
        #
        npcmti[isp] = lmmaxi*nrcmti[isp]
        npcmt[isp] = npcmti[isp] + lmmaxo*(nrcmt[isp] - nrcmti[isp])
    end


    SHT = SphericalHarmonicTransform(lmaxi, lmaxo)

    gntyry = _init_gntyry(lmaxapw, idxlm, lmaxo)

    return MuffinTins(
        nrmtscf, nrmt, rmtall, rmtdelta, rmt, omegamt, lradstp,
        nrcmt, rcmt, rlmt, rlcmt, wrmt, wprmt, wrcmt, wprcmt,
        maxlapw, lmaxapw, lmmaxapw, lmaxo, lmmaxo, lmaxi, lmmaxi, fracinr,
        nrmti, nrcmti, idxlm, idxil, idxim, npmti, npmt, npcmti, npcmt,
        SHT, gntyry
    )

end

function _init_gntyry(lmaxapw, idxlm, lmaxo)
    lmmaxo = (lmaxo+1)^2
    lmmaxapw = (lmaxapw+1)^2
    gntyry = zeros(ComplexF64, lmmaxo, lmmaxapw, lmmaxapw)
    for l1 in 0:lmaxapw, m1 in -l1:l1
        lm1 = idxlm[l1,m1]
        for l3 in 0:lmaxapw, m3 in -l3:l3
            lm3 = idxlm[l3,m3]
            for l2 in 0:lmaxo, m2 in -l2:l2
                lm2 = idxlm[l2,m2]
                gntyry[lm2,lm3,lm1] = gauntyry(l1,l2,l3,m1,m2,m3)
            end
        end
    end
    return gntyry
end


# TO BE REMOVED
function MuffinTins(Nspecies::Int64; lmaxi=1)

    nrmtscf = 0
    nrmt = zeros(Int64,Nspecies)
    rmtall = 0
    rmtdelta = 0.05 # default
    rmt = zeros(Nspecies)
    omegamt = 0.0
    lradstp = 4
    nrcmt = zeros(Int64,Nspecies)

    # XXX SHould be done in genrmesh
    rcmt = Vector{Vector{Float64}}(undef,Nspecies)
    rlmt = Vector{OffsetArray{Float64,2,Array{Float64,2}}}(undef,Nspecies)
    rlcmt = Vector{OffsetArray{Float64,2,Array{Float64,2}}}(undef,Nspecies)
    
    wrmt = Vector{Vector{Float64}}(undef,Nspecies)
    wprmt = Vector{Matrix{Float64}}(undef,Nspecies)
    wrcmt = Vector{Vector{Float64}}(undef,Nspecies)
    wprcmt = Vector{Matrix{Float64}}(undef,Nspecies)
    
    maxlapw = 50
    lmaxapw = 8 # default
    lmmaxapw = (lmaxapw+1)^2
    lmaxo  = 6
    lmmaxo = (lmaxo+1)^2
    lmaxi  = min(lmaxi,lmaxo)
    lmmaxi = (lmaxi+1)^2
    fracinr = 0.01
    
    nrmti = zeros(Int64,Nspecies)
    nrcmti = zeros(Int64,Nspecies)

    # index to (l,m) pairs
    idxlm = OffsetArray( zeros(Int64,lmaxapw+1,2*lmaxapw+1),
        0:lmaxapw,-lmaxapw:lmaxapw)
    idxil = zeros(Int64, lmmaxapw)
    idxim = zeros(Int64, lmmaxapw)
    lm = 0
    for l in 0:lmaxapw, m in -l:l
        lm = lm + 1
        idxlm[l,m] = lm
        idxil[lm] = l
        idxim[lm] = m
    end

    npmti = zeros(Int64,Nspecies)
    npmt = zeros(Int64,Nspecies)
    npcmti = zeros(Int64,Nspecies)
    npcmt = zeros(Int64,Nspecies)

    SHT = SphericalHarmonicTransform(lmaxi, lmaxo)

    return MuffinTins(
        nrmtscf, nrmt, rmtall, rmtdelta, rmt, omegamt, lradstp,
        nrcmt, rcmt, rlmt, rlcmt, wrmt, wprmt, wrcmt, wprcmt,
        maxlapw, lmaxapw, lmmaxapw, lmaxo, lmmaxo, lmaxi, lmmaxi, fracinr,
        nrmti, nrcmti, idxlm, idxil, idxim, npmti, npmt, npcmti, npcmt,
        SHT
    )

end

# TO BE REMOVED
# why need this separate step?
function init_zero!( mt_vars::MuffinTins )
    
    nrmt = mt_vars.nrmt
    Nspecies = size(nrmt, 1)

    nrcmt = mt_vars.nrcmt
    lradstp = mt_vars.lradstp
    
    # make the muffin-tin mesh commensurate with lradstp
    for is in 1:Nspecies
        nrmt[is] = nrmt[is] - (nrmt[is]-1)%lradstp
        nrcmt[is] = (nrmt[is] - 1)/lradstp + 1
    end

    @assert mt_vars.lmaxo <= mt_vars.lmaxapw

    return
end


# Should be called after genrmesh
# This initializes the following fields:
#   npmti, npmt
#   npcmti, npcmt
function init_packed_mtr!( mt_vars::MuffinTins )
    #
    nrmti = mt_vars.nrmti
    Nspecies = size(nrmti,1)
    nrmt = mt_vars.nrmt
    #
    nrcmt = mt_vars.nrcmt
    nrcmti = mt_vars.nrcmti
    #
    lmmaxo = mt_vars.lmmaxo
    lmmaxi = mt_vars.lmmaxi
    #

    for isp in 1:Nspecies
        #
        mt_vars.npmti[isp] = lmmaxi*nrmti[isp]
        mt_vars.npmt[isp] = mt_vars.npmti[isp] + lmmaxo*(nrmt[isp] - nrmti[isp])
        #
        mt_vars.npcmti[isp] = lmmaxi*nrcmti[isp]
        mt_vars.npcmt[isp] = mt_vars.npcmti[isp] + lmmaxo*(nrcmt[isp] - nrcmti[isp])
    end
    return
end