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

# Using default value for xmin, rmax, zmesh
function RadialGrid( Nrmesh::Int64 )
    r = zeros(Float64, Nrmesh)
    r2 = zeros(Float64, Nrmesh)
    rab = zeros(Float64, Nrmesh)
    sqrtr = zeros(Float64, Nrmesh)
    rm1 = zeros(Float64, Nrmesh)
    rm2 = zeros(Float64, Nrmesh)
    rm3 = zeros(Float64, Nrmesh)
    xmin = 0.0
    rmax = 0.0
    zmesh = 0.0
    dx = 0.0
    return RadialGrid(Nrmesh, r, r2, rab, sqrtr, rm1, rm2, rm3, xmin, rmax, zmesh, dx)
end


import Base: show
function show(io::IO, pw::RadialGrid)
    @printf(io, "RadialGrid object\n")
    return
end

# Build the radial (logarithmic) grid 
#
# r(i) = exp ( xmin + (i-1) dx ) / zmesh  i=1,mesh
# r2(i) is r(i) square, sqr(i) is sqrt(r(i)) and 
# rab(i) is the integration element = r(i)*dx
#
# more general grid definitions are possible but currently not implemented
# (example: Vanderbilt's grid, same as above but starting at r=0)
# r(i) = exp ( xmin ) * ( exp( (i-1)*dx ) - 1.0_dp ) / zmesh
# rab(i) = ( r(i) + exp(xmin)/zmesh ) * dx


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
    # mesh must be odd for Simpson integration
    Nrmesh = round(Int64, 2*floor(Nrmesh/2) + 1 )

    if Nrmesh > NrmeshMax
        println("WARNING: Nrmesh is larger than NrmeshMax")
    end
    
    if ibound
        xmin = xmax - dx*(Nrmesh-1)
    end

    r = zeros(Float64, Nrmesh)
    r2 = zeros(Float64, Nrmesh)
    rab = zeros(Float64, Nrmesh)
    sqrtr = zeros(Float64, Nrmesh)
    rm1 = zeros(Float64, Nrmesh)
    rm2 = zeros(Float64, Nrmesh)
    rm3 = zeros(Float64, Nrmesh)

    for ir in 1:Nrmesh
        x = xmin + (ir-1)*dx # linear grid
        r[ir] = exp(x)/zmesh
        r2[ir] = r[ir]*r[ir]
        rab[ir] = r[ir]*dx
        sqrtr[ir] = sqrt(r[ir])
        rm1[ir] = 1.0/r[ir]
        rm2[ir] = 1.0/r[ir]^2
        rm3[ir] = 1.0/r[ir]^3
    end
    
    return RadialGrid(Nrmesh, r, r2, rab, sqrtr, rm1, rm2, rm3, xmin, rmax, zmesh, dx)
end

# simple routine returning the coefficient of the polynomial 
# describing the leading behavior of a function f at small r.
#
# TODO: some bibliography
function radial_grid_series!(f, r, r2, b)


    dr21 = r[2] - r[1]
    dr31 = r[3] - r[1]
    dr32 = r[3] - r[2]
    dr41 = r[4] - r[1]
    dr42 = r[4] - r[2]
    dr43 = r[4] - r[3]
    df21 = (f[2] - f[1])/dr21
    df32 = (f[3] - f[2])/dr32
    df43 = (f[4] - f[3])/dr43
    ddf42 = (df43 - df32)/dr42
    ddf31 = (df32 - df21)/dr31

    # Note that b is originally indexed as 0:3
    # we convert it to 1-based indexing
    b[4] = (ddf42 - ddf31)/dr41
    b[3] = ddf31 - b[4]*( r[1] + r[2] + r[3] )
    b[2] = df21 - b[3]*( r[2] + r[1] ) - b[4]*( r2[1] + r2[2] + r[1]*r[2] )
    b[1] = f[1] - r[1]*( b[2] + r[1]*(b[3] + r[1]*b[4]) )

   return
end
