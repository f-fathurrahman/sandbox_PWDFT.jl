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
  return RadialGrid(Nrmesh, r2, rab, sqrtr, rm1, rm2, rm3, xmin, rmax, zmesh, dx)
end

