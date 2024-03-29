from ase.units import kJ, Bohr, Hartree
from ase.eos import EquationOfState
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 150

filename = "TEMP_EOS_data.dat"
dat = np.loadtxt(filename)

# dont forget to convert to Angstrom and eV
volumes = dat[:,0]**3/4.0 * Bohr**3  # only for fcc
energies = dat[:,1]*Hartree
energies_abinit = dat[:,2]*Hartree
energies_pwscf  = dat[:,3]*Hartree

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print("PWDFT.jl: B = %18.10f GPa" % (B/kJ * 1.0e24))
eos.plot("pwdft-eos-ase.png")

eos = EquationOfState(volumes, energies_abinit)
v0, e0, B = eos.fit()
print("ABINIT  : B = %18.10f GPa" % (B/kJ * 1.0e24))
eos.plot("abinit-eos-ase.png")

eos = EquationOfState(volumes, energies_pwscf)
v0, e0, B = eos.fit()
print("PWSCF   : B = %18.10f GPa" % (B/kJ * 1.0e24))
eos.plot("pwscf-eos-ase.png")
