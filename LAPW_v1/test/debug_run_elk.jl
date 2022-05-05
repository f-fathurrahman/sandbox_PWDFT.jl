using OffsetArrays
using Serialization: serialize

# import PyPlot as plt

const LIBLAPW = "/home/efefer/WORKS/my_github_repos/ffr-PWDFT/src_LAPW/liblapwdft.so"

ccall( (:read_input_, LIBLAPW), Cvoid, () )
ccall( (:init0_, LIBLAPW), Cvoid, () )
ccall( (:init1_, LIBLAPW), Cvoid, () )
ccall( (:info_gvectors_, LIBLAPW), Cvoid, () )
ccall( (:info_muffin_tins_, LIBLAPW), Cvoid, () )
ccall( (:writesym_, LIBLAPW), Cvoid, () )
ccall( (:my_rhoinit_, LIBLAPW), Cvoid, () )
ccall( (:maginit_, LIBLAPW), Cvoid, () )
ccall( (:my_potks_, LIBLAPW), Cvoid, (Ref{Bool},), true )
