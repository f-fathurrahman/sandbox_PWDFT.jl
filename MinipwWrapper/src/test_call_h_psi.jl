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
psi = zeros(ComplexF64, npwx, nbnd)
for i in 1:nbnd
    psi[i,i] = 1.0
end
hpsi = zeros(ComplexF64, npwx, nbnd)
spsi = zeros(ComplexF64, npwx, nbnd)

# XXX: Note that this will allocate some unmanaged memory
# All references to becmod (and related modules and types) should be modified
# in the Fortran subroutines
ccall( (:prepare_h_s_psi_, LIBMINIPW), Cvoid, (Ref{Int32},), Int32(ik) )


println("Before: sum(psi) = ", sum(psi))
println("Before: sum(hpsi) = ", sum(hpsi))
println("Before: sum(spsi) = ", sum(spsi))

# SUBROUTINE my_h_psi( lda, n, m, psi, hpsi )
ccall(
    (:my_h_psi_, LIBMINIPW), Cvoid, 
    (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ptr{ComplexF64}),
    Int32(npwx), Int32(ngk[ik]), Int32(nbnd), psi, hpsi
)

# SUBROUTINE my_s_psi( lda, n, m, psi, spsi )
ccall(
    (:my_s_psi_, LIBMINIPW), Cvoid, 
    (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ptr{ComplexF64}),
    Int32(npwx), Int32(ngk[ik]), Int32(nbnd), psi, spsi
)

println("After: sum(psi) = ", sum(psi))
println("After: sum(hpsi) = ", sum(hpsi))
ss = ( sum(abs.(spsi)) - sum(psi) )*0.25
println("After: sum(spsi) = ", ss) # to Ha (?)
