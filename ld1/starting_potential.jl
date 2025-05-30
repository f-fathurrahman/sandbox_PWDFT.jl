#
# Copyright (C) 2004 PWSCF group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt.
#
# starting potential: Generalized Thomas-Fermi atomic potential
# it is assumed that the effective charge seen by the reference
# electron cannot be smaller than 1 (far from the core)
# TODO: add reference
#
function starting_potential!(
  Nrmesh, Zval, Zed, Nwf, oc, nn, ll, r, enl,
  v0, vxt, Vpot;
  frozen_core=false, noscf=false
)
  
    enne = 0.0
    zz = max(Zed, Zval)
    
    for n in 1:Nwf
       oce = max(0.0, oc[n])
       enne = enne + oce
       zen= 0.0
       for i in 1:Nwf
            oce = max(0.0, oc[i])
            if nn[i] < nn[n] 
                zen = zen + oce
            end
            if (nn[i] == nn[n]) && (ll[i] <= ll[n])
                zen = zen + oce
            end
        end
        zen = max(zz - zen + 1.0, 1.0)
        if (abs(enl[n]) < 1.e-7) || (!frozen_core)
            enl[n] = -0.5*( zen/nn[n] )^2  # Ha unit
        end
    end
    
    for i in 1:Nrmesh
       #vxt[i] = Vext( r[i] )
       vxt[i] = 0.0
       x = r[i]*enne^(1.0/3.0)/0.885
       t = zz/(1.0 + sqrt(x)*(0.02747 - x*(0.1486 - 0.007298*x)) + x*(1.243 + x*(0.2302 + 0.006944*x)))
       t = max(1.0,t)
       v0[i] = -Zed/r[i]
       if noscf
          Vpot[i] = v0[i] + vxt[i]
       else
          Vpot[i] = -t/r[i] + vxt[i]
       end
    end
    
    # XXX: spinpol is handled outside this function
    return
end
