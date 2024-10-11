
#=
Returns the Fermi-Dirac approximation to the Dirac delta function
$$ \tilde\delta(x)=\frac{e^{-x}}{(1+e^{-x})^2}. $$
=#
function sdelta_fd(x)
    if abs(x) > 50.0
        return 0.0
    end
    t1 = exp(-x)
    return t1/( (1.0 + t1)^2 )
end

