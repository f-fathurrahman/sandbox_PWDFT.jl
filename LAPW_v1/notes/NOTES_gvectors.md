## Generation of G-vectors

In Elk, to generate G-vectors we must specify gmaxvr.

We will use ecutrho instead for this purpose. Note that:

```julia
ecutrho = 0.5*gmaxvr^2
```

`gmaxvr` must be at least `2*gkmax`

In Elk, the variable `gc` saves the norm of G-vectors (not the norm-squared).

`ngridg` is the variable that corresponds to Ns in PWDFT.jl (grid size
along the 1st, 2nd, and 3rd lattice vectors).

ngtot is total grid FFT points (product of ngridg)

## Generation of Gk-vectors

Elk also use gkmax instead of ecutwfc for initializing Gk-vectors.

`gkmax` is determined by quantity rgkmax, which is the muffin tin radius times gkmax. Elk can use largest, smallest, average or user specified muffin tin radius. The determination of gkmax is done in before generation of G-vectors in `init_gvector_arrays`.

## Atomic structure factor

This is done in subroutine gensfacgp().

The atomic structure factors for a set of ${\bf G+p}$-vectors are
defined by:
$$
S_{\alpha}({\bf G+p})=\exp(i({\bf G+p})\cdot{\bf r}_{\alpha})
$$
where ${\bf r}_{\alpha}$ is the position of atom $\alpha$. 

This is different convention from what we used in PWDFT.jl where we used negative sign in front of imaginary number.

## Smooth characteristic function

Done by subroutine gencfun(). The variable is cfunig(ngtot) and cfunir(ngtot).

cfunig will be used in olpistl.f90

cfunir will be used in moment.f90.

Generates the smooth characteristic function. This is the function which is 0 within the muffin-tins and 1 in the intersitial region and is constructed from radial step function form factors with $G<G_{{\rm max}}$.


