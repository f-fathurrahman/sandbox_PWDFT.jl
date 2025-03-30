# Tobe copy-pasted to REPL for various debugging stuffs

using DelimitedFiles

# read formatted Fortran files
d_r = readdlm("compare_qe/Si_v01/fort.1111");
Npoints = length(d_r);
d_r = reshape(d_r, (Npoints,));

e_r = readdlm("compare_qe/Si_v01/fort.1112");
Npoints = length(e_r);
e_r = reshape(e_r, (Npoints,));

rhs_r = readdlm("compare_qe/Si_v01/fort.1113");
Npoints = length(rhs_r);
rhs_r = reshape(rhs_r, (Npoints,));

V_h_r = readdlm("compare_qe/Si_v01/fort.1114");
Npoints = length(V_h_r);
V_h_r = reshape(V_h_r, (Npoints,));

