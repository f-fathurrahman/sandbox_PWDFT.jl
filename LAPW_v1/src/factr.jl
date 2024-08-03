function factr(n::Int64, d::Int64)
    @assert n >= 0
    @assert d >= 0
    if d == 1
        return factnm(n,1)
    end
    if n < d
        res = Float64(n+1)
        for i in (n+2):d
            res = res*i
        end
        return 1.0/res
    elseif n == d
        return 1.0
    else
        res = Float64(d+1)
        for i in (d+2):n
            res = res*i
        end
        return res
    end
end




