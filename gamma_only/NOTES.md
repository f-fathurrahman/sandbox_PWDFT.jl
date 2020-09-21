# Gamma-point Trick Implementation Notes

## Overview

Gamma-point trick is used for real quantities that depend of $\mathbf{G}$-vectors.
Using the trick, only (about) half of the $\mathbf{G}$-vectors are stored.
The expansion coefficient for other half can be calculated from:
$$
c(\mathbf{G}) = c^{*}(-\mathbf{G})
$$

The $\mathbf{G}$-vectors are now described with `GVectorsGamma` instead of `GVectors`.
From `GVectorsGamma` we can derive `GVectorsWGamma`.
Here is the definition of `GVectorsGamma`:

```julia
struct GVectorsGamma
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
    idx_g2rm::Array{Int64,1}
    G2_shells::Array{Float64,1}
    idx_g2shells::Array{Int64,1}
end
```

Only about half of the actual $\mathbf{G}$-vectors components are stored.
There is now a new field, namely `idx_g2rm` which stores the indices of 'negative'
of $\mathbf{G}$-vectors.

## Generating GVectorsGamma

The algorithm for generating half set of $\mathbf{G}$-vectors are adapted from the
`ggen` subroutine of PWSCF. I simplify a bit and use adapt it according to the convention
that I have used before.

Mapping from 3d index to linear index.

## Operations involving wavefunction

Notes: I have not made any comparison with PWSCF implementation. I come up with my own
implementation by studying some operations with special arrays.

Test array:

dot product

overlap matrix

Wavefunction orthonormalization: ortho_GS_gamma and ortho_sqrt_gamma.


FIXME: Need to compare gradient calculation: calc_grad using Gamma-trick vs usual.


## Calculating electron density

Blah:

## Hamiltonian

New struct instead of parametric  struct: `HamiltonianGamma`. Several things that change

Construction of local potential