# tolerance for error in total charge
function get_epschg()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_epschg, LIBLAPW), Float64))
end

# total nuclear charge
function get_chgzn()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgzn, LIBLAPW), Float64))
end

# total valence charge
function get_chgval()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgval, LIBLAPW), Float64))
end

# total charge
function get_chgtot()
    return unsafe_load(cglobal((:__m_charge_moment_current_MOD_chgtot, LIBLAPW), Float64))
end


#=
# core charges
REAL(8) chgcr(maxspecies)
# total core charge
REAL(8) chgcrtot
# core leakage charge
REAL(8), ALLOCATABLE :: chgcrlk(:)
# excess charge
REAL(8) chgexs
# calculated total charge
REAL(8) chgcalc
# interstitial region charge
REAL(8) chgir
# muffin-tin charges
REAL(8), ALLOCATABLE :: chgmt(:)
# total muffin-tin charge
REAL(8) chgmttot
# effective Wigner radius
REAL(8) rwigner
# total moment
REAL(8) momtot(3)
# total moment magnitude
REAL(8) momtotm
# interstitial region moment
REAL(8) momir(3)
# muffin-tin moments
REAL(8), ALLOCATABLE :: mommt(:,:)
# total muffin-tin moment
REAL(8) mommttot(3)
# total current
REAL(8) curtot(3)
# total current magnitude
REAL(8) curtotm
=#