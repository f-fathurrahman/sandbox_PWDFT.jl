function r3frac!(v; epslat=1e-6)
    for i in 1:3
        v[i] = v[i] - round(Int64, v[i])
        if v[i] < 0.0
            v[i] = v[i] + 1.0
        end
        if (1.0 - v[i]) < epslat
            v[i] = 0.0
        end
        if v[i] < epslat
            v[i] = 0.0
        end
    end
    return
end 