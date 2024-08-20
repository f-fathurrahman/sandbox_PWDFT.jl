const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )
#ccall( (:my_electrons_, LIBMINIPW), Cvoid, () )

npwx = unsafe_load(cglobal((:__wvfct_MOD_npwx, LIBMINIPW), Int32)) |> Int64
nbnd = unsafe_load(cglobal((:__wvfct_MOD_nbnd, LIBMINIPW), Int32)) |> Int64
npol = unsafe_load(cglobal((:__noncollin_module_MOD_npol, LIBMINIPW), Int32)) |> Int64

nks = unsafe_load(cglobal((:__klist_MOD_nks, LIBMINIPW), Int32)) |> Int64
println("nks = ", nks)
# Read ngk
ptr = cglobal( (:__klist_MOD_ngk, LIBMINIPW), Ptr{Int32} )
ngk = zeros(Int64, nks)
for i in 1:nks
    ngk[i] = unsafe_load(unsafe_load(ptr,1),i) |> Int64
end

println("npwx = ", npwx)
println("ngk = ", ngk)
println("nbnd = ", nbnd)


ik = 1

# XXX: Note that this will allocate some unmanaged memory
# All references to becmod (and related modules and types) should be modified
# in the Fortran subroutines
ccall( (:prepare_h_s_psi_, LIBMINIPW), Cvoid, (Ref{Int32},), Int32(ik) )

Ngw_k = ngk[ik]
H_diag = zeros(Float64, Ngw_k)
S_diag = zeros(Float64, Ngw_k)

# SUBROUTINE my_usnldiag(npw, h_diag, s_diag)
ccall( (:my_usnldiag_, LIBMINIPW), Cvoid,
    (Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
    Int32(Ngw_k), H_diag, S_diag
)

println("sum H_diag = ", sum(H_diag)*0.5) # to Ha
println("sum S_diag = ", sum(S_diag) - Ngw_k)
