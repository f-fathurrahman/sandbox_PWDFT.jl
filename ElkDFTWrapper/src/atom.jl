function get_natmtot()
    natmtot = unsafe_load(cglobal( (:__m_atoms_MOD_natmtot, LIBLAPW), Int32 )) |> Int64
    return natmtot
end

function elk_solve_atom!(
    SOL, ptnucl,
    zn, nst_,
    n_, l_, k_, occ,
    xctype_, xcgrad,
    nr_, r,
    evals, rho, vr, rwf
)

# Arguments with trailing underscore must be converted to Int32

#  REAL(8), intent(in) :: sol
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