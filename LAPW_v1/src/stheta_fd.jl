#=
Returns the Fermi-Dirac approximation to the Heaviside step function
$$ \tilde\Theta(x)=\frac{1}{1+e^{-x}}. $$
=#
function stheta_fd(x::Float64)
    # large
    if x > 50.0
        return 1.0
    end
    # small
    if x < -50.0
        return 0.0
    end
    # default
    return 1.0/(1.0 + exp(-x))
end
