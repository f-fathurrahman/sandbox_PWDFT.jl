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
end

function MuffinTins(Nspecies; lmaxi=1)

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
    for l in 0:lmaxapw
        for m in -l:l
            lm = lm + 1
            idxlm[l,m] = lm
            idxil[lm] = l
            idxim[lm] = m
        end
    end

    npmti = zeros(Int64,Nspecies)
    npmt = zeros(Int64,Nspecies)
    npcmti = zeros(Int64,Nspecies)
    npcmt = zeros(Int64,Nspecies)

    return MuffinTins(
        nrmtscf, nrmt, rmtall, rmtdelta, rmt, omegamt, lradstp,
        nrcmt, rcmt, rlmt, rlcmt, wrmt, wprmt, wrcmt, wprcmt,
        maxlapw, lmaxapw, lmmaxapw, lmaxo, lmmaxo, lmaxi, lmmaxi, fracinr,
        nrmti, nrcmti, idxlm, idxil, idxim, npmti, npmt, npcmti, npcmt
    )

end


function init_zero!( mt_vars::MuffinTins )
    
    nrmt = mt_vars.nrmt
    Nspecies = size(nrmt)[1]

    nrcmt = mt_vars.nrcmt
    lradstp = mt_vars.lradstp
    
    # make the muffin-tin mesh commensurate with lradstp
    for is in 1:Nspecies
        nrmt[is] = nrmt[is] - (nrmt[is]-1)%lradstp
        nrcmt[is] =( nrmt[is] - 1)/lradstp + 1
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