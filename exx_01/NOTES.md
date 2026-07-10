## NOTES

In `pw.x`, EXX can be activated by using `input_dft='hf-noc'` or some other
functionals such as `pbe0` or `hse`.

In `qe-6.6`, using `input_dft='hf-nox-noc'` is still OK, but `qe-7.4` (and some
previous versions) will complain.

I found out that, in the first inner (?) iteration (when EXX is not active),
input_dft='hf-noc' will actually use Slater exchange, which solely depends
on density, please see subroutine `xc` in `xc_lda_lsda_drivers.f90` in QE.

In order to completely disabled XC (pure HF via EXX?), I need to modify the
subroutine such that it returns zeros.

There is some differences when I start with xc none and Slater xc. This
is happened when I use smearing, the one starting Slater xc (the original QE code)
is giving correct (?) total energy, with smearing contribution zero. Meanwhile
starting from xc none give non zero smearing contribution.

## Some steps

(1) Try to build exx grid (q-grid)
