using Infiltrator

using SpecialFunctions
using LinearAlgebra
using Dates
using Random

import Optim
import LineSearches


includet("MV.jl")
includet("wannierize_utils.jl")
includet("run_optim.jl")

function debug_run_optim()

    filename = "TEMP_free"
    read_amn = true #read $file.amn as input
    read_eig = true #read the eig file (can be set to false when not disentangling)
    do_write_amn = true #write $file.optimize.amn at the end
    nfrozen = 1 #will freeze if either n <= nfrozen or eigenvalue in frozen window
    frozen_window_low = -Inf
    frozen_window_high = -Inf
    ftol = 1e-20 #tolerance on spread
    gtol = 1e-4 #tolerance on gradient
    maxiter = 3000 #maximum optimization iterations
    m = 100 #history size of BFGS

    # expert/experimental features
    do_normalize_phase = false # perform a global rotation by a phase factor at the end
    do_randomize_gauge = false #randomize initial gauge
    cluster_size = 1e-6 #will also freeze additional eigenvalues if the freezing cuts a cluster. Set to 0 to disable
    only_r2 = false #only minimize sum_n <r^2>_n, not sum_n <r^2>_n - <r>_n^2

    Random.seed!(0)

    p = read_system(filename, read_amn, read_eig)


    if read_amn
        A0 = copy(p.A)
    else
        A0 = randn(size(p.A)) + im*randn(size(p.A))
    end

    for i in 1:p.N1, j in 1:p.N2, k in 1:p.N3
        l_frozen, l_not_frozen = local_frozen_sets(
            p, nfrozen, i, j, k,
            frozen_window_low, frozen_window_high,
            cluster_size
        )
        A0[:,:,i,j,k] = normalize_and_freeze(A0[:,:,i,j,k], l_frozen, l_not_frozen)
    end

    A = minimize(
        p, A0,
        nfrozen, frozen_window_low, frozen_window_high,
        cluster_size,
        do_randomize_gauge
    )

    # fix global phase
    if do_normalize_phase
        for i in 1:nwannier
            imax = indmax(abs.(A[:,i,1,1,1]))
            @assert abs(A[imax,i,1,1,1]) > 1e-2
            A[:,i,:,:,:] *= conj(A[imax,i,1,1,1]/abs(A[imax,i,1,1,1]))
        end
    end

    if do_write_amn
        write_amn(p,A,"$(p.filename).optimized")
    end


    loc = omega_loc(p,A)
    # writedlm("data",loc)
    println(maximum(loc))
    println("Max loc = ", maximum(loc))
    #include("plot_free_wannier.jl")
end