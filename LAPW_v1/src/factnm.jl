# Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
# This file is distributed under the terms of the GNU Lesser General Public
# License. See the file COPYING for license details.
# !REVISION HISTORY:
#   Created January 2003 (JKD)

# !INPUT/OUTPUT PARAMETERS:
#   n : input (in,integer)
#   m : order of multifactorial (in,integer)
# !DESCRIPTION:
#   Returns the multifactorial
#   $$ n\underbrace{!!\,\cdots\,!}_{m\,{\rm times}}=
#    \prod_{\substack{i\ge 0\\ n-im>0}}(n-im) $$
#   for $n,\,m \ge 0$. $n$ should be less than 150.

function factnm(n, m)
    f1 = [
                           1.0,                        2.0,
                           6.0,                       24.0,
                         120.0,                      720.0,
                        5040.0,                    40320.0,
                      362880.0,                  3628800.0,
                    39916800.0,                479001600.0,
                  6227020800.0,              87178291200.0,
               1307674368000.0,           20922789888000.0,
             355687428096000.0,         6402373705728000.0,
          121645100408832000.0,      2432902008176640000.0,
        51090942171709440000.0,   1124000727777607680000.0,
     25852016738884976640000.0, 620448401733239439360000.0]

    f2 = [
                           1.0,                        2.0,
                           3.0,                        8.0,
                          15.0,                       48.0,
                         105.0,                      384.0,
                         945.0,                     3840.0,
                       10395.0,                    46080.0,
                      135135.0,                   645120.0,
                     2027025.0,                 10321920.0,
                    34459425.0,                185794560.0,
                   654729075.0,               3715891200.0,
                 13749310575.0,              81749606400.0,
                316234143225.0,            1961990553600.0,
               7905853580625.0,           51011754393600.0,
             213458046676875.0,         1428329123020800.0,
            6190283353629375.0,        42849873690624000.0,
          191898783962510625.0,      1371195958099968000.0,
         6332659870762850625.0,     46620662575398912000.0,
       221643095476699771875.0,   1678343852714360832000.0,
      8200794532637891559375.0,  63777066403145711616000.0]
    
    # fast return if possible
    if n == 0
        return 1.0
    end

    if m == 1
        if (n >= 1) && (n <= 24)
          return f1[n]
        end
    end

    if m == 2
        if (n >= 1) && (n <= 38)
            return f2[n]
        end
    end

    if n < 0
      error("Error(factnm): n < 0")
    end

    if m <= 0
        error("Error(factnm): m <= 0")
    end
    
    if n > 150
        error("Error(factnm): n out of range")
    end

    if m == 1
        res = f1[24]
        for i in 25:n
          res = res*i
        end
        return res
    else
        j = floor(Int64, n/m)
        if mod(n,m) == 0
            j = j - 1
        end
        res = Float64(n)
        for i in 1:j
            res = res*(n-i*m)
        end
        return res
    end

    return 0.0 # FIXME to avoid type-unstability?
end

# Compare with fortran:
# ccall( (:factnm_, LIBLAPW), Float64, (Ref{Int32}, Ref{Int32}), Int32(20), Int32(3))