using PWDFT

pw = PWGrid(15.0, gen_lattice_sc(5.0))
println(pw)

G = pw.gvec.G
Ng = pw.gvec.Ng
lmax = 2
lmmax = (lmax + 1)^2
dylm = zeros(Float64, Ng, lmmax)
dYlm_real_qe!(lmax, G, dylm, 1)

