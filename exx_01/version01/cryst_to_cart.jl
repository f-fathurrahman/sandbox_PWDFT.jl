#ffr: This is not intuitive
#ffr: It seems better to use RecVecs and LatVecs with inv although it might be slower
#
# will modify v
function cryst_to_cart!(nv, v, trmat, iflag)
    @assert iflag in [-1,1]
    vau = zeros(Float64, 3)
    for nv in 1:nv
        if iflag == 1
            for k in 1:3
                vau[k] = trmat[k,1]*v[1,nv] + trmat[k,2]*v[2,nv] + trmat[k,3]*v[3,nv]
            end
        else
            for k in 1:3
                vau[k] = trmat[1,k]*v[1,nv] + trmat[2,k]*v[2,nv] + trmat[3,k]*v[3,nv]
            end
        end
        for k in 1:3
            v[k,nv] = vau[k]
        end
    end
    return
end