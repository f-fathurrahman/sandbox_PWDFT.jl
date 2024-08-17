# ElkDFTWrapper

A simple wrapper to Elk DFT package.

Mainly used for writing variables or arrays to Julia's serialized data.

Some convention (hopefully I consistently used this):

- Elk's original subroutines are all prefixed with function `call_` followed with
  the subroutine name. This is similar to subroutine call syntax in Fortran.
  Previously, `elk_` prefix is used, however, the package name already contains Elk
  string, so this is somewhat redundant.

- "Getter" functions are not exported.
  They must be called with `ElkDFTWrapper.getter_function`.

- No specific convention for other function names.

## Hamiltonian and overlap matrices construction

Inputs:
- radial basis functions
- plane waves
- potentials (muffin tins and interstitial)

Relevant subroutines:
- hmlfv, olpfv
- Called by hmlfv:
  - hmlaa, olpaa (for each atoms)
  - hmlistl, olpistl (for interstitital)
  - hmlalo, olpalo (for all atoms)
  - hmllolo, olplolo (for all atoms)
