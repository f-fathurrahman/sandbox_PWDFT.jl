using ElkDFTWrapper
const elk = ElkDFTWrapper

elk.init_debug_calc()
vkc = elk.get_vkc()
ik = 3
@info "vkc = $(vkc[:,ik])"
apwalm, Hmat, Omat = elk.get_hamiltonian(ik=ik, ispin=1)

