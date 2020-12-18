mutable struct RadialGrid
    Nrmesh::Int64           # the actual number of mesh points
    r::Array{Float64,1}     # the radial mesh
    r2::Array{Float64,1}    # the square of the radial mesh
    rab::Array{Float64,1}   # d r(x) / d x where x is the linear grid
    sqrtr::Array{Float64,1} # the square root of the radial mesh
    rm1::Array{Float64,1}   # 1 / r
    rm2::Array{Float64,1}   # 1 / r^2
    rm3::Array{Float64,1}   # 1 / r^3
    xmin::Float64       # the minimum x
    rmax::Float64       # the maximum radial point
    zmesh::Float64      # the ionic charge used for the mesh
    dx::Float64         # the deltax of the linear mesh
end

function RadialGrid( Nrmesh::Int64 )
    r = zeros(Float64,Nrmesh)
    r2 = zeros(Float64,Nrmesh)
    rab = zeros(Float64,Nrmesh)
    sqrtr = zeros(Float64,Nrmesh)
    rm1 = zeros(Float64,Nrmesh)
    rm2 = zeros(Float64,Nrmesh)
    rm3 = zeros(Float64,Nrmesh)
    xmin = 0.0
    rmax = 0.0
    zmesh = 0.0
    dx = 0.0
    return RadialGrid(Nrmesh, r, r2, rab, sqrtr, rm1, rm2, rm3, xmin, rmax, zmesh, dx)
end


# from do_mesh of upflib
function RadialGrid(
    rmax::Float64,
    zmesh::Float64,
    xmin::Float64,
    dx::Float64,
    ibound::Bool
)

    NrmeshMax = 3500 # parameter from radial_grid, probably not needed
    
    xmax = log(rmax*zmesh)
    
    Nrmesh = round(Int64, floor((xmax-xmin)/dx) + 1)
    # mesh must be odd for simpson integration
    Nrmesh = round(Int64, 2*floor(Nrmesh/2) + 1 )

    if Nrmesh > NrmeshMax
        println("WARNING: Nrmesh is larger than NrmeshMax")
    end
    
    if ibound
        xmin = xmax - dx*(Nrmesh-1)
    end

    r = zeros(Float64,Nrmesh)
    r2 = zeros(Float64,Nrmesh)
    rab = zeros(Float64,Nrmesh)
    sqrtr = zeros(Float64,Nrmesh)
    rm1 = zeros(Float64,Nrmesh)
    rm2 = zeros(Float64,Nrmesh)
    rm3 = zeros(Float64,Nrmesh)

    for ir in 1:Nrmesh
        x = xmin + (ir-1)*dx
        r[ir]   = exp(x)/zmesh
        r2[ir]  = r[ir]*r[ir]
        rab[ir] = r[ir]*dx
        sqrtr[ir] = sqrt(r[ir])
        rm1[ir] = 1.0/r[ir]
        rm2[ir] = 1.0/r[ir]^2
        rm3[ir] = 1.0/r[ir]^3
    end
    
    return RadialGrid(Nrmesh, r, r2, rab, sqrtr, rm1, rm2, rm3, xmin, rmax, zmesh, dx)
end