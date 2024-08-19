using ElkDFTWrapper
const elk = ElkDFTWrapper

elk.init_debug_calc()

vxcmt = elk.get_vxcmt();
vxcmt_orig = vxcmt[:,:];

npmtmax = elk.get_npmtmax()
nrmt = elk.get_nrmt()
nrmti = elk.get_nrmti()
symlatc = elk.get_symlatc()

isp = 1
ia = 1
rot = symlatc[:,:,2]
@views rfmt1 = vxcmt[:,ia];
rfmt2 = zeros(Float64, npmtmax);
elk.call_rotrfmt!(rot, nrmt[isp], nrmti[isp], rfmt1, rfmt2)


