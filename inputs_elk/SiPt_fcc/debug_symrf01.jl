using ElkDFTWrapper
const elk = ElkDFTWrapper

elk.init_debug_calc()

vxcmt = elk.get_vxcmt()
vxcir = elk.get_vxcir()

vxcmt_orig = vxcmt[:,:]
vxcir_orig = vxcir[:]

npmtmax = elk.get_npmtmax()
elk.call_symrf!(vxcmt, vxcir)


