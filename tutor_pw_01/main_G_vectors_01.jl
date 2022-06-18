using LinearAlgebra: diagm
using PWDFT

# Copy of init_gvec
function init_idx_miller( Ns, RecVecs, ecutrho )

    Ng = PWDFT.calc_Ng( Ns, RecVecs, ecutrho )

    G_temp = zeros(Float64,3)
    G2 = Array{Float64}(undef,Ng)
    idx_miller = Array{Int64}(undef,3,Ng)

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1, j in 0:Ns[2]-1, i in 0:Ns[1]-1
        ip = ip + 1
        gi = PWDFT.mm_to_nn( i, Ns[1] )
        gj = PWDFT.mm_to_nn( j, Ns[2] )
        gk = PWDFT.mm_to_nn( k, Ns[3] )
        G_temp[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G_temp[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G_temp[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2_temp = G_temp[1]^2 + G_temp[2]^2 + G_temp[3]^2
        if 0.5*G2_temp <= ecutrho
            ig = ig + 1
            G2[ig] = G2_temp
            idx_miller[1,ig] = gi
            idx_miller[2,ig] = gj
            idx_miller[3,ig] = gk
        end
    end

    idx_sorted = sortperm(G2)
    idx_miller = idx_miller[:,idx_sorted]

    return idx_miller
end


LatVecs = gen_lattice_fcc(3.0)
#LatVecs = gen_lattice_tetragonal_P(4.0, 2.0)
display(LatVecs); println()

h = LatVecs
a1 = LatVecs[:,1]
a2 = LatVecs[:,2]
a3 = LatVecs[:,3]

display(a1); println()
display(a2); println()
display(a3); println()

# Calculate unit reciprocal vectors: b1, b2, b3
RecVecs = Matrix(2π*inv(h'))
display(RecVecs); println()
println(typeof(RecVecs))

ecutwfc = 2.0
pw = PWGrid(ecutwfc, LatVecs)
println(pw)

idx_miller = init_idx_miller(pw.Ns, pw.RecVecs, pw.ecutrho)
println(size(idx_miller))

for ig in 1:pw.gvec.Ng
    Gvec = RecVecs*idx_miller[:,ig]
    #println("Gvec = ", Gvec, " pw.gvec.G = ", pw.gvec.G[:,ig])
    println("ΔGvec = ", Gvec - pw.gvec.G[:,ig])
end
