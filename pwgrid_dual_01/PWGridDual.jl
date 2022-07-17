const PLANFW_TYPE = PWDFT.PLANFW_TYPE
const PLANBW_TYPE = PWDFT.PLANBW_TYPE

struct PWGridDual
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    Nss::Union{Tuple{Int64,Int64,Int64},Nothing}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    CellVolume::Float64
    gvec::GVectors
    gvecs::Union{GVectors,Nothing}
    gvecw::GVectorsW
    planfw::PLANFW_TYPE
    planbw::PLANBW_TYPE
    planfws::Union{PLANFW_TYPE,Nothing}
    planbws::Union{PLANBW_TYPE,Nothing}
end


function PWGridDual(
    ecutwfc::Float64, LatVecs::Array{Float64,2};
    kpoints=nothing, Ns_=(0,0,0), dual=4.0
)

    @assert dual >= 4.0

    ecutrho = dual*ecutwfc
    RecVecs = 2*pi*inv(Matrix(LatVecs'))
    CellVolume = abs(det(LatVecs))
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])


    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1
    if any(Ns_ .== 0)
        Ns1 = PWDFT.good_fft_order(Ns1)
        Ns2 = PWDFT.good_fft_order(Ns2)
        Ns3 = PWDFT.good_fft_order(Ns3)
        Ns = (Ns1,Ns2,Ns3)
    else
        Ns = Ns_[:]
    end
    Npoints = prod(Ns)
    gvec = PWDFT.init_gvec( Ns, RecVecs, ecutrho )
    #
    planfw = plan_fft!( zeros(ComplexF64,Ns) )
    planbw = plan_ifft!( zeros(ComplexF64,Ns) )

    if dual > 4.0
        Ns1 = 2*round( Int64, sqrt(4*ecutwfc/2)*LatVecsLen[1]/pi ) + 1
        Ns2 = 2*round( Int64, sqrt(4*ecutwfc/2)*LatVecsLen[2]/pi ) + 1
        Ns3 = 2*round( Int64, sqrt(4*ecutwfc/2)*LatVecsLen[3]/pi ) + 1
        Ns1 = PWDFT.good_fft_order(Ns1)
        Ns2 = PWDFT.good_fft_order(Ns2)
        Ns3 = PWDFT.good_fft_order(Ns3)
        Nss = (Ns1,Ns2,Ns3)
        # XXX Setting Nss from the constructor is not yet supported
        #
        # TODO: simply copy from gvec instead of calling this function again
        # Need to calculate Ngs (no. of G-vectors for smooth grid)
        gvecs = PWDFT.init_gvec( Nss, RecVecs, 4*ecutwfc )
        planfws = plan_fft!( zeros(ComplexF64,Nss) )
        planbws = plan_ifft!( zeros(ComplexF64,Nss) )
    else
        Nss = nothing
        gvecs = nothing
        planfws = nothing
        planbws = nothing
    end

    if kpoints == nothing
        kpoints = KPoints( 1, (1,1,1), zeros(3,1), [1.0], RecVecs )
    end

    gvecw = PWDFT.init_gvecw( ecutwfc, gvec, kpoints )

    return PWGridDual(
        ecutwfc, ecutrho, Ns, Nss, LatVecs, RecVecs, CellVolume,
        gvec, gvecs, gvecw, planfw, planbw, planfws, planbws
    )

end


import Base: show
function show( io::IO, pw::PWGridDual; header=true )
    if header
        @printf("\n")
        @printf("                                 ----------\n")
        @printf("                                 PWGridDual\n")
        @printf("                                 ----------\n")
        @printf("\n")
    end
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    @printf(io, "Direct lattice vectors:\n")
    @printf(io, "\n")
    for i = 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf(io, "\n")
    @printf(io, "Reciprocal lattice vectors:\n")
    @printf(io, "\n")
    for i = 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf(io, "\n")
    @printf(io, "Direct lattive volume = %18.10f bohr^3\n", pw.CellVolume )
    @printf(io, "ecutwfc               = %18.10f Ha\n", pw.ecutwfc)
    @printf(io, "ecutrho               = %18.10f Ha\n", pw.ecutrho)    
    @printf(io, "Sampling points       = (%5d,%5d,%5d)\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    #
    show( io, pw.gvec )
    show( io, pw.gvec, pw.gvecw )
end
show( pw::PWGrid; header=true ) = show( stdout, pw, header=header )