# Read CIF
import ase.io
atoms = ase.io.read("TaAs_mp-1936_primitive.cif")

# Reduce the cell, also modify the atomic coordinates
from ase.build import niggli_reduce
niggli_reduce(atoms)

# Export coordinate to xyz and POSCAR
atoms.write("TaAs_primitive.xyz")
atoms.write("POSCAR_TaAs", format="vasp")

from ase.dft.kpoints import get_special_points
spec_kpts = get_special_points(atoms.cell)
print(spec_kpts)
