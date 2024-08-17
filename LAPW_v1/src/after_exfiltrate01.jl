# After calling @exfiltrate
# To be pasted into REPL


# For testing calc_match_coeffs
atoms = exfiltrated.atoms;
mt_vars = exfiltrated.mt_vars;
apwlo_vars = exfiltrated.apwlo_vars;
pw = exfiltrated.pw;
ik = 1
apwalm = Vector{Array{ComplexF64,3}}(undef, atoms.Natoms);
for ia in 1:atoms.Natoms
    isp = atoms.atm2species[ia]
    apwordmax = maximum(apwlo_vars.apword[isp])
    Ngk = pw.gvecw.Ngw[ik]
    lmmaxapw = mt_vars.lmmaxapw
    apwalm[ia] = zeros(ComplexF64, Ngk, apwordmax, lmmaxapw);
end

calc_match_coeffs!(ik, atoms, pw, mt_vars, apwlo_vars, apwalm);

# For testing Hamiltonian construction
#
haa = exfiltrated.haa;
nmat = exfiltrated.nmat;
H = zeros(ComplexF64, nmat[ik], nmat[ik]);
for ia in 1:atoms.Natoms
    hmlaa!(ik, ia, atoms, pw, mt_vars, apwlo_vars, apwalm, haa, H);
end