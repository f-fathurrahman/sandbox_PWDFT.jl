function call_rhoinit()
    ccall( (:rhoinit_, LIBLAPW), Cvoid, () )
    return
end

function call_potks(; txc=true)
    ccall( (:potks_, LIBLAPW), Cvoid, (Ref{Bool},), true )
    return
end

function call_potcoul()
    ccall( (:potcoul_, LIBLAPW), Cvoid, () )
    return
end

# We call potxc_default here instead of potxc as it is simpler
function call_potxc()
    ccall( (:potxc_default_, LIBLAPW), Cvoid, () )
    return
end


# TODO: make a better interface for this
function call_atom!(
    SOL, ptnucl,
    zn, nst_,
    n_, l_, k_, occ,
    xctype_, xcgrad,
    nr_, r,
    evals, rho, vr, rwf
)
# SOL: 

#  REAL(8), intent(in) :: sol speed of light (in atomic unit)
#  LOGICAL, intent(in) :: ptnucl
#  REAL(8), intent(in) :: zn
#  INTEGER, intent(in) :: nst
#  INTEGER, intent(in) :: n(nst),l(nst),k(nst)
#  REAL(8), intent(inout) :: occ(nst)
#  INTEGER, intent(in) :: xctype(3),xcgrad
#  INTEGER, intent(in) :: nr
#  REAL(8), intent(in) :: r(nr)
#  REAL(8), intent(out) :: eval(nst)
#  REAL(8), intent(out) :: rho(nr),vr(nr)
#  REAL(8), intent(out) :: rwf(nr,2,nst)

    # NOTE: Arguments with trailing underscore must be converted to Int32
    nst = convert(Int32, nst_)
    n = convert(Array{Int32}, n_)
    l = convert(Array{Int32}, l_)
    k = convert(Array{Int32}, k_)
    nr = convert(Int32, nr_)
    xctype = convert(Array{Int32}, xctype_)

    ccall(
        (:atom_, LIBLAPW),
        Cvoid, (
            Ref{Float64}, Ref{Bool},
            Ref{Float64}, Ref{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
            Ptr{Int32}, Ref{Bool},
            Ref{Int32}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        SOL, ptnucl,
        zn, nst,
        n, l, k, occ,
        xctype, xcgrad,
        nr, r,
        evals, rho, vr, rwf
    )

    return
end
