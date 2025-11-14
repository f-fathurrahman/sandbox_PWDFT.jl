using Infiltrator, Serialization

atoms = exfiltrated.atoms;
mt_vars = exfiltrated.mt_vars;
apwlo_vars = exfiltrated.apwlo_vars;
pw = exfiltrated.pw;
elec_chgst = exfiltrated.elec_chgst;

nstfv = elec_chgst.nstfv;
ik = 1;
ispin = 1;
#evecfv = deserialize("evecs_1st_ispin_$(ispin)_ik_$(ik).dat");
#apwalm = calc_match_coeffs(ik, atoms, pw, mt_vars, apwlo_vars);

evecfv = deserialize("evecs_1st_ispin_$(ispin)_ik_$(ik)_elk.dat");
apwalm_ia = deserialize("apwalm_ispin_1_ik_1_elk.dat")

ia = 1; isp = 1;
npcmtmax = maximum(mt_vars.npcmt);
lradstp = mt_vars.lradstp;

wfmt1 = zeros(ComplexF64, npcmtmax, nstfv);

ist = 1;
fill!(wfmt1, 0.0);
println("sum evecfv[:,ist] = ", sum(evecfv[:,ist]))
@views wavefmt!( lradstp, ia, atoms, mt_vars, apwlo_vars,
    pw.gvecw.Ngw[ik], apwalm_ia, evecfv[:,ist], wfmt1[:,ist] )
println("sum wfmt1[:,ist] = ", sum(wfmt1[:,ist]))


# compute the first-variational wavefunctions
for ist in 1:nstfv
    println("sum evecfv[:,ist] = ", sum(evecfv[:,ist]))
    @views wavefmt!( lradstp, ia, atoms, mt_vars, apwlo_vars,
        pw.gvecw.Ngw[ik], apwalm_ia, evecfv[:,ist], wfmt1[:,ist] )
    println("sum wfmt1[:,ist] = ", sum(wfmt1[:,ist]))
end


ndmag = exfiltrated.ndmag; # should be 1?
bsmt = exfiltrated.bsmt;
nsc = 2;
wrcmt = mt_vars.wrcmt;

wfmt2 = zeros(ComplexF64, npcmtmax);
wfmt3 = zeros(ComplexF64, npcmtmax);
wfmt4 = zeros(ComplexF64, npcmtmax, nsc);
npc = mt_vars.npcmt[isp];
println("sum(wfmt1) = ", sum(wfmt1));
println("sum bsmt[ia][1:npc,ndmag] = ", sum(bsmt[ia][1:npc,ndmag]));

jst = 1;
# convert wavefunction to spherical coordinates
fill!(wfmt2, 0.0);
println("sum(wfmt1[:,jst]) = ", sum(wfmt1[:,jst]));
@views backward_SHT!(mt_vars, isp, wfmt1[:,jst], wfmt2, coarse=true);
println("sum(wfmt2) = ", sum(wfmt2))
#
# apply Kohn-Sham effective magnetic field
fill!(wfmt3, 0.0);
wfmt3[1:npc] .= bsmt[ia][1:npc,ndmag] .* wfmt2[1:npc];  # ffr: THIS IS IMPORTANT
println("sum(wfmt3) = ", sum(wfmt3))
#
# convert to spherical harmonics and store in wfmt4
fill!(wfmt4, 0.0);
@views forward_SHT!(mt_vars, isp, wfmt3, wfmt4[:,1], coarse=true);
wfmt4[1:npc,2] .= -wfmt4[1:npc,1]
println("sum(wfmt4[:,1]) = ", sum(wfmt4[:,1]))
println("sum(wfmt4[:,2]) = ", sum(wfmt4[:,2]))

#
z1 = zfmtinp(mt_vars, isp, wrcmt[isp], wfmt1[:,ist], wfmt4, coarse=true)
println("Upper diagonal block: ist=$ist z1 from zfmtinp = $z1")

