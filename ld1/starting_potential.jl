#
# Copyright (C) 2004 PWSCF group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
#

function starting_potential!(
  Vext,
  Nrmesh, Zval, Zed, Nwf, oc, nn, ll, r, enl,
  v0, vxt, vpot, enne, nspin;
  frozen_core=false, noscf=false
)
    #
    # starting potential: generalized thomas-fermi atomic potential
    # it is assumed that the effective charge seen by the reference
    # electron cannot be smaller than 1 (far from the core)
    #
    #implicit none
    #integer :: nwf, nn(nwf), ll(nwf), ndm, mesh, n, i, nspin
    #real(DP) :: r(ndm), vpot(ndm,2), v0(ndm), vxt(ndm), enl(nwf), oc(nwf), &
    #   zed, zval, zz, zen, enne, t,x, vext, oce
    #real(DP), parameter :: e2 = 2.0
    #external vext
  
    enne = 0.0
    zz = max(Zed, Zval)
    
    for n in 1:Nwf
       oce = max(0.0, oc[n])
       enne = enne + oce
       zen= 0.0
       for i = 1:nwf
            oce = max(0.0, oc[i])
            if nn[i] < nn[n] 
                zen = zen + oce
            end
            if (nn[i] == nn[n]) && (ll[i] <= ll[n])
                zen = zen + oce
            end
        end
        zen = max( zz - zen + 1.0, 1.0 )
        if ABS( enl(n)) < 1.e-7 || (!frozen_core)
            enl[n] = -( zen/nn[n] )^2
        end
    end
    
    for i in 1:Nrmesh
       vxt[i] = Vext( r[i] )
       x = r[i]*enne^(1.0/3.0)/0.885
       t = zz/(1.0 + sqrt(x)*(0.02747 - x*(0.1486 - 0.007298*x)) + x*(1.243 + x*(0.2302 + 0.006944*x)))
       t = max(1.0,t)
       v0[i]= -zed/r[i]
       if noscf
          vpot[i,1] = v0[i] + vxt[i]
       else
          vpot[i,1] = -t/r[i] + vxt[i]
       end
    end
    
    if Nspin == 2
       for i = 1:Nrmesh
          vpot[i,2] = vpot[i,1]
       end
    end
    return
end
