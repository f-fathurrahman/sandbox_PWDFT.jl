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
psi = ones(ComplexF64, npwx, nbnd)
hpsi = zeros(ComplexF64, npwx, nbnd)
spsi = zeros(ComplexF64, npwx, nbnd)

# XXX: Note that this will allocate some unmanaged memory
# All references to becmod (and related modules and types) should be modified
# in the Fortran subroutines
ccall( (:prepare_h_s_psi_, LIBMINIPW), Cvoid, (Ref{Int32},), Int32(ik) )


println("Calling h_psi")
ccall(
    (:my_h_psi_, LIBMINIPW), Cvoid, 
    (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{ComplexF64}, Ptr{ComplexF64}),
    Int32(npwx), Int32(ngk[ik]), Int32(nbnd), psi, hpsi
)
println("End of calling h_psi")
# FIXME: call to calbec should be sufficient

ccall( (:__jl_comm_MOD_copy_becp_k, LIBMINIPW), Cvoid, (), )

ptr = cglobal( (:__jl_comm_MOD_becp_k, LIBMINIPW), Ptr{ComplexF64} )
Ndim1 = unsafe_load(cglobal((:__jl_comm_MOD_nrow_becp_k, LIBMINIPW), Int32)) |> Int64
Ndim2 = unsafe_load(cglobal((:__jl_comm_MOD_ncol_becp_k, LIBMINIPW), Int32)) |> Int64
tmp = zeros(ComplexF64, Ndim1*Ndim2)
ip = 1
for j in 1:Ndim2, i in 1:Ndim1
    tmp[ip] = unsafe_load(unsafe_load(ptr,1),ip)
    ip = ip + 1
end
becp_k = reshape(tmp, Ndim1, Ndim2)
println("sum becp_k = ", sum(becp_k)*0.5)

display(abs.(becp_k[1:5,1:5])*0.5); println();
