const LIBMINIPW = "/home/efefer/WORKS/my_github_repos/ffr-learns-qe/minipw-6.6/srclib/libqemain.so"

ccall( (:prepare_all_, LIBMINIPW), Cvoid, () )
ccall( (:info_upf_, LIBMINIPW), Cvoid, () )
ccall( (:my_electrons_, LIBMINIPW), Cvoid, () )
#ccall( (:my_forces_, LIBMINIPW), Cvoid, () )
#stress = zeros(Float64,3,3)
#ccall( (:my_stress_, LIBMINIPW), Cvoid, (Ptr{Float64},), stress )
#println("\nStress obtained from my_stress\n")
#display(stress); println()
