# !INPUT/OUTPUT PARAMETERS:
#   n  : number of required points (in,integer)
#   tp : (theta, phi) coordinates (out,real(2,n))
function sphcover!(n::Int64, tp::Array{Float64,2})
    @assert n > 0
    dz = 2.0/n
    z = 1.0 - dz/2.0
    tp[1,1] = acos(z)
    dp = pi*(1.0 - sqrt(5.0))
    p = 0.0
    tp[2,1] = p
    for k in 2:n
        z = z - dz
        tp[1,k] = acos(z)
        p = p + dp
        # rem here does the same thing as mod in f90
        # Use rem instead of mod
        tp[2,k] = rem(p, 2*pi)
    end
    return
end 
