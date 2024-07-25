# No spherical restriction
# G2 shells are not implemented
struct GVectorsFull
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
    idx_g2miller::Vector{Tuple{Int64,Int64,Int64}}
end

function GVectorsFull( Ns::Tuple{Int64,Int64,Int64}, RecVecs::Matrix{Float64} )

    Ng = prod(Ns) # only for full gvectors

    G_temp = zeros(Float64,3)

    G  = Array{Float64}(undef,3,Ng)
    G2 = Array{Float64}(undef,Ng)
    idx_g2r = Array{Int64}(undef,Ng)
    idx_g2miller = Vector{Tuple{Int64,Int64,Int64}}(undef,Ng) # Array{Int64,2}(undef,3,Ng)

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1, j in 0:Ns[2]-1, i in 0:Ns[1]-1
        ip = ip + 1
        gi = PWDFT.mm_to_nn( i, Ns[1] )
        gj = PWDFT.mm_to_nn( j, Ns[2] )
        gk = PWDFT.mm_to_nn( k, Ns[3] )
        # XXX RecVecs are stored by rows
        G_temp[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G_temp[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G_temp[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2_temp = G_temp[1]^2 + G_temp[2]^2 + G_temp[3]^2
        # No check here if G2_temp is less than some value
        ig = ig + 1
        @views  G[:,ig] = G_temp[:]
        G2[ig] = G2_temp
        idx_g2r[ig] = ip
        idx_g2miller[ig] = (gi, gj, gk)
    end

    # Sort G vector according to its magnitude
    idx_sorted = sortperm(G2)
    G = G[:,idx_sorted]
    G2 = G2[idx_sorted]
    idx_g2r = idx_g2r[idx_sorted]
    idx_g2miller = idx_g2miller[idx_sorted]
    return GVectorsFull( Ng, G, G2, idx_g2r, idx_g2miller )

end
