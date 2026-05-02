 Jerzy Goraus - 2011-05-18

Dear elk users and developers,
I have a question regarding the basis, let's take for example 'Mn.in'

```text
'Mn'                                       : spsymb
'manganese'                                : spname
  -25.0000                                  : spzn
   100145.9369                              : spmass
  0.400000E-06    2.0000   33.4391   500    : sprmin, rmt, sprmax, nrmt
  10                                        : spnst
   1   0   1   2.00000    T                 : spn, spl, spk, spocc, spcore
   2   0   1   2.00000    T
   2   1   1   2.00000    T
   2   1   2   4.00000    T
   3   0   1   2.00000    F
   3   1   1   2.00000    F
   3   1   2   4.00000    F
   3   2   2   3.00000    F
   3   2   3   2.00000    F
   4   0   1   2.00000    F
   1     : apword
  0.1500   0  F  : apwe0, apwdm, apwve
   0    : nlx
   6    : nlorb
   0   2                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
   1   2                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
   2   2                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
   0   3                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
 -3.3872   0  T
   1   3                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
 -2.2357   0  T
   2   3                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
 -0.5052   0  T
```

My question is which states belong to APW+lo basis and which to LAPW+LO.
If I understand correctly apword 1  means that we have somehow APW basis (this energy goes into  formula 5.1
in S. Cottenier  introduction to (L)APW from 2004) .
Then we have zero exceptions from above the APW basis,
and then we have six local orbitals.

What concerns the three first local orbitals they look like formula 5.2 for APW+lo or 4.2 from LAPW, this 0.15 Ha
is the E1,l linearization energy.  But then we have three another local orbitals which
look like LAPW+LO (formula 4.4) as they have E1,l  and E2,l and also derivative term.
If they would be APW+lo+LO then they would lack the derivative term (formula 5.3).
But as stated above by the nlx  we have no exceptions from APW+lo ?

So my qustions are:
- are the first three states LAPW or APW+lo ?
- are the last three states LAPW+LO ?
- in the valence section we have that 3s, 3p,3d and 4s are valence/semicore states, how to understand
that  (LAPW or APW+lo) and LAPW+LO is defined for the same states, as the lines with the same lorbl are
repeated ?

sincerely,
Jerzy Goraus 

---

 Lars Nordström - 2011-05-18

Dear Jerzy,

OK I agree, the nomenclature is a bit confusing and the input file of ELK has a larger flexibility than what is mirrored by the conventional acronyms.

The APW+lo and LAPW+LO refer to that there are essentially two types of basis functions.
So the first thing to understand is the difference of an augmented plane wave (APW/LAPW)and a local orbital (lo/LO)
They simply differ in the interstitial region, the APW is a plane wave and the LO vanishes.

However inside a muffin-tin both these are expanded in local spherical functions in essentially the same way.
For each angular momentum at each atom, the radial part can be expressed with one or several radial functions. Each radial functions corresponds to a partial wave solution at a certain energy or to an energy derivate of a certain order of that function at a given energy.

Now the lorbord parameter is the number of different radial functions for the given angular momentum lorbl.
However, in principle you can mix the energies and derivatives in any way you want for each lo (or apw). Hence the forms you find in Cottenier's description are only special cases, though the two most common.

Now you do not want to use many radial functions in the APW part, most often only one, due to slower convergence.
When you use two in the APW, for historically reasons it is referred to as a LAPW …

In the example you have given. The basis consists of:
1) APW's with one radial function (apword=1) for each atom and each l, with no exception (nlx=0). The number of APW's is determined by the PW-cutoff of the parameter, rgkmax.
and
2) six independent lo's
The first three with two radial functions. A function and its derivative at same energy, in all three cases.
The three last lo's use three radial functions. A function and its derivative at same energy plus another function at a lower energy.

The first three lo's are to cover the valence states, ie 4s, 4p and 3d.
the fourth and fifth cover the semicores, 3s, 3p,
while the last one add variational freedom for the 3d states in order to treat the fairly large exchange splitting of Mn.

In Wien2k-terminology (as used by Cottenier) the six lo's would best be described as: 3 lo, 2 LO and 1 lo, respectively.
However, the 6'th lo has three radial functions and therefore is in form close to the LO, but is not intended for semi-core states … I prefer to use "lo" for all kinds of lo's in order to avoid such confusion.

This input file is by no means unique, you can alter it in many ways and get the same result. If one wants to understand the basis set better that is the way to go! I give three examples below how to get your fingers dirty:

For instance I suggest that you try to use two exception (nlx=2), with two radial functions for the 4s and 4p partial waves and then skip the corresponding lo's, ending up with four lo's. This results in a somewhat smaller total basis set but probably give very similar results. Why? This would end up in a mixed APW/LAPW basis set. Therefore I prefer to use APW for all possible choices, other than the classical LAPW choice, with two radial function, the partial wave and its energy derivative at one given linearisation energy.

You never really have to go higher order (lorbord) than two. How would that look like in the Mn.in file?

Another thing is that in principle we never need to use energy derivatives, but instead partial waves at different energies. Try to rewrite the Mn.in file in order to do that. The energy derivatives are there for completeness and historical reasons, as they were introduced in the linearisation of the old APW-method, i.e. in the LAPW.

Best regards,
    Lars Nordström 

