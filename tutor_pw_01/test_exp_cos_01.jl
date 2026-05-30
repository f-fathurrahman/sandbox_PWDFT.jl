function exp_cos(x)
    α = 0.5
    x0 = 0.5
    T = 2.0
    arg1 = (2π/T) * (x - x0)
    return exp(-α) * exp(α*cos(arg1))
end

function gaussian_cos(x)
    α = 0.5
    x0 = 0.5
    T = 2.0
    arg1 = (2π/T) * (x - x0)
    return exp(-α^2) * exp(α^2*cos(arg1))
end

