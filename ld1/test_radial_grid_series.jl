include("RadialGrid.jl")

function main()
    f = [1.0, 2.0, 3.0, 4.0]
    r = [1.1, 1.2, 1.3, 1.4]
    r2 = r.^2
    b = zeros(Float64,4)
    radial_grid_series!(f, r, r2, b)
    println("b = ", b)
end

main()