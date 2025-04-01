# FIXME: probably can be replaced by LAPACK
function solve_tridiag!(a, b, c, r, u, n)

    if abs(b[1]) < 1.e-10
        error("b(1) is too small")
    end

    γ = zeros(Float64, n)
    β = b[1]
    u[1] = r[1]/β
    for j in 2:n
        γ(j) = c[j-1]/β
        β = b[j] - a[j] * γ[j]
        if abs(β) < 1.e-10
            error("β is too small")
        end
        u[j] = ( r[j] - a[j]*u[j-1] )/β
    end
    for j in range(n-1,stop=1,step=-1)
        u[j] = u[j] - γ[j+1] * u[j+1]
    end
    return
end