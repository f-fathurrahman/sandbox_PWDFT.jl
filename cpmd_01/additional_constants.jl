# From QE
# Should be moved to PWDFT

const AMU_SI = 1.660538782e-27
const ELECTRONMASS_SI = 9.10938215e-31
const AMU_AU = AMU_SI / ELECTRONMASS_SI

const H_PLANCK_SI = 6.62607015e-34      # J s
const HARTREE_SI = 4.3597447222071e-18  # J
const AU_SEC = H_PLANCK_SI/(2*pi)/HARTREE_SI
const AU_PS  = AU_SEC*1.0e12