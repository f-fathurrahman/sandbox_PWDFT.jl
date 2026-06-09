using Infiltrator

using SpecialFunctions
using LinearAlgebra
using Dates

includet("write_files.jl")

#=
We look at the free Hamiltonian H = -Delta in the unit cell 2pi x [-.5,.5]^d, with BZ [0,1]^d
For every k, we therefore solve the Hamiltonian H_k (-i nabla + k)^2 on 2pi [-.5,.5]^d with periodic BC
This is discretized in the basis of K/G vectors e_K(r) = exp(iK.r), with K integers
In this basis, H_KK' = delta_KK' |K|^2
Therefore unk(r) = e^(i K_n r), where K_n corresponds to the n-th element
of the K ordered by increasing value
if perm = sortperm(eig_k), then perm gives the mapping from n to K
Ki is limited to -Li:Li.
By Fourier duality this gives a real-space grid of the unit cell

i,j,k are indices in the BZ (from 1 to N1,N2,N3), kind is the flat version of (i,j,k) (from 1 to N1*N2*N3).
The BZ loops are always on [0,1]^d, which is shifted by "corner" when
computing λ, so that the actual BZ is corner + [0,1]^d
=#

# all loops in the code refer to a k point grid like (0:N-1)/N.
# This transfers to the actual kpt grid
function k_red_to_real(k, NN; corner = [0.0, 0.0, 0.0], do_shift = false)
    N1, N2, N3 = NN # unpack
    if do_shift
        shift = [1/N1/2, 1/N2/2, 1/N3/2]
        shift[1:3-dim] = 0 #don't shift singleton dimensions
        return k + shift + corner
    else
        return k + corner
    end
end

function k_real_to_red(k, NN; corner = [0.0, 0.0, 0.0], do_shift = false)
    N1, N2, N3 = NN # unpack
    if do_shift
        shift = [1/N1/2, 1/N2/2, 1/N3/2]
        shift[1:3-dim] = 0 # don't shift singleton dimensions
        return k - shift - corner
    else
        return k - corner
    end
end

# overlap matrix Mmnkb between k and k+b, with k+b = displ_vec + kpb
# <un(k)|um(kpb + displ_vec)> = <un(k)|e^-(displ_vec)x um(kpb)>
# = delta_KK', but mn are sorted versions of K K'
function get_ovl(k, kpb, displ_vec, Ks, NN, nband)
    eig_k = get_eigvals(k, Ks, NN)
    eig_kpb = get_eigvals(kpb, Ks, NN)
    sort_k = sortperm(eig_k)
    sort_kpb = sortperm(eig_kpb)
    A = zeros(Int64, size(Ks,2),size(Ks,2))
    for m in 1:nband, n in 1:nband
        Km = sort_k[m]
        @assert eig_k[Km] == sort(eig_k)[m]
        Kpbn = sort_kpb[n]
        A[m,n] = Ks[:,Km] == (Ks[:,Kpbn] .- displ_vec)
    end
    return A
end

function get_eigvals(k, Ks, NN)
    # eigenvalues at k in [0,1]^d
    return dropdims(sum(abs2, (Ks .+ k_red_to_real(k, NN)), dims=1), dims=1)
end


function debug_main()

    nwannier = 2
    # k-grid
    NN = [1, 1, 20]
    N1, N2, N3 = NN

    # G-grid
    LL = [0, 0, 2]
    L1, L2, L3 = LL

    dim = 3 - count(x -> x==1, NN)
    println("dim = $dim")

    #SCDM parameters
    mu = 0.0
    sigma = 2.0

    do_shift = false
    corner = [0.,0.,0.]
    # corner = [0.,0.,0.] # do not use directly in the code, only the functions below

    k = randn(3);
    @assert k ≈ k_red_to_real(k_real_to_red(k, NN), NN) #sanity check

    nband = (2L1+1)*(2L2+1)*(2L3+1) #size of the Hamiltonian
    Ks = zeros(Int, 3, nband)
    let ind = 1
        for l1 in -L1:L1, l2=-L2:L2, l3 in -L3:L3
            Ks[:,ind] = [l1, l2, l3]
            ind += 1
        end
    end
    println("nband = $nband")

    write_files(NN, LL, nwannier, nband, mu, sigma)

    return

end