## Elk Species file

This file is needed for initialization of atomic species, MuffinTins
and APLWO.

Probably it is better to use separate files for reading, for example:
Cu.species, Cu.radial_mt, Cu.apwlo

For now, it might be better to use readspecies as is. Need to read
genspecies to tweak the configuration.

Example of species file:

```
 'Cu'                                       : spsymb
 'copper'                                   : spname
  -29.0000                                  : spzn
   115837.2717                              : spmass
  0.371391E-06    2.4000   42.9740   500    : rminsp, rmt, rmaxsp, nrmt
  10                                        : nstsp
   1   0   1   2.00000    T                 : nsp, lsp, ksp, occsp, spcore
   2   0   1   2.00000    T
   2   1   1   2.00000    T
   2   1   2   4.00000    T
   3   0   1   2.00000    T
   3   1   1   2.00000    F
   3   1   2   4.00000    F
   3   2   2   4.00000    F
   3   2   3   6.00000    F
   4   0   1   1.00000    F
   1                                        : apword
    0.1500   0  F                           : apwe0, apwdm, apwve
   0                                        : nlx
   5                                        : nlorb
   0   2                                    : lorbl, lorbord
    0.1500   0  F                           : lorbe0, lorbdm, lorbve
    0.1500   1  F
   1   2                                    : lorbl, lorbord
    0.1500   0  F                           : lorbe0, lorbdm, lorbve
    0.1500   1  F
   2   2                                    : lorbl, lorbord
    0.1500   0  F                           : lorbe0, lorbdm, lorbve
    0.1500   1  F
   1   2                                    : lorbl, lorbord
    0.1500   0  F                           : lorbe0, lorbdm, lorbve
   -2.6152   0  T
   2   2                                    : lorbl, lorbord
    0.1500   0  F                           : lorbe0, lorbdm, lorbve
   -0.1918   0  T
```

The variable nlx seems to be always 0.

