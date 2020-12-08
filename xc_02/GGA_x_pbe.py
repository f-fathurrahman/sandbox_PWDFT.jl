from sympy import *
from numeric_symbols import *
#from symbolic_symbols import *

X2S        = 1/(2*(6*pi**2)**(1/3))

X_FACTOR_C = THREE/EIGHT*(THREE/pi)**(ONE/THREE) * 4**(TWO/THREE)
# Dimension = 3
RS_FACTOR = (THREE/(FOUR*pi))**(ONE/THREE)
LDA_X_FACTOR = -X_FACTOR_C
DIMENSIONS = 3

def lda_x_spin(rs, z):
    return LDA_X_FACTOR*((1 + z)/2)**(1 + 1/DIMENSIONS)*(RS_FACTOR/rs)

# Polarization = ferr
def gga_exchange(func, rs, z, xs0, xs1):
    return lda_x_spin(rs, 1)*func(xs0)

#def gga_exchange(func, rs, z, xs0, xs1):
#    return lda_x_spin(rs, z)*func(xs0) + lda_x_spin(rs, -z)*func(xs1)

# (* standard PBE *)
#$ifdef gga_x_pbe_params
params_a_kappa = 0.8040
params_a_mu    = 0.2195149727645171
#$endif

#(* PBE_SOL *)
#$ifdef gga_x_pbe_sol_params
#params_a_kappa := 0.8040:
#params_a_mu    := MU_GE:
#$endif

#$ifdef gga_x_pbe_tca_params
#params_a_kappa := 1.227:
#params_a_mu    := 0.2195149727645171:
#$endif

def pbe_f0(s):
    return 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + params_a_mu*s**2))

def pbe_f(x):
    return pbe_f0(X2S*x)

def f(rs, z, xt, xs0, xs1):
    return gga_exchange(pbe_f, rs, z, xs0, xs1)

r_s, z, xt, xs0, xs1 = symbols("r_s z xt xs0 xs1")
ε_x = f(r_s, z, xt, xs0, xs1)

from sympy.utilities.codegen import codegen

code1 = codegen( ("eps_x", ε_x), language="julia")
print(code1[0][1])

d_ε_x1 = diff(ε_x, r_s)
code1 = codegen( ("d_eps_x1", d_ε_x1), language="julia")
print(code1[0][1])

d_ε_x2 = diff(ε_x, xs0)
code1 = codegen( ("d_eps_x2", d_ε_x2), language="julia")
print(code1[0][1])
