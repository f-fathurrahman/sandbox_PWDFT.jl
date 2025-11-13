using Serialization

atoms = exfiltrated.atoms;
mt_vars = exfiltrated.mt_vars;
apwlo_vars = exfiltrated.apwlo_vars;
pw = exfiltrated.pw;
elec_chgst = exfiltrated.elec_chgst;

nstfv = elec_chgst.nstfv;
ik = 1;
ispin = 1;
evecfv = deserialize("evecs_1st_ispin_$(ispin)_ik_$(ik).dat");
apwalm = calc_match_coeffs(ik, atoms, pw, mt_vars, apwlo_vars);

ia = 1;
npcmtmax = maximum(mt_vars.npcmt);
wfmt1 = zeros(ComplexF64, npcmtmax, nstfv);

lradstp = mt_vars.lradstp;
ist = 1;
@views wavefmt!( lradstp, ia, atoms, mt_vars, apwlo_vars,
    pw.gvecw.Ngw[ik], apwalm[ia], evecfv[:,ist], wfmt1[:,ist] )

# compute the first-variational wavefunctions
for ist in 1:nstfv
    @views wavefmt!( lradstp, ia, atoms, mt_vars, apwlo_vars,
        pw.gvecw.Ngw[ik], apwalm[ia], evecfv[:,ist], wfmt1[:,ist] )
end