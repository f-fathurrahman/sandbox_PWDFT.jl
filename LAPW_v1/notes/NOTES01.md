Several implemented data structures:
```
src/APWLOVars.jl:mutable struct APWLOVars
src/AtomicSpeciesVars.jl:mutable struct AtomicSpeciesVars
src/AtomicVars.jl:struct AtomicVars
src/LatticeVars.jl:struct LatticeVars
src/Libxc_old.jl:mutable struct XCFuncType
src/MuffinTins.jl:mutable struct MuffinTins
src/SphericalHarmonicTransform.jl:struct SphericalHarmonicTransform
src/SymmetryVars.jl:mutable struct SymmetryVars
src/XCCalculator.jl:struct XCCalculator <: AbstractXCCalculator
src/XCCalculator.jl:struct LibxcXCCalculator <: AbstractXCCalculator
```

They are converted from global variables in Elk. Fortunately, Elk has
grouped these global variables, so we converted them to `struct`s.

