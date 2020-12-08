from sympy import *
#from numeric_symbols import *
from symbolic_symbols import *

X_FACTOR_C = THREE/EIGHT*(THREE/pi)**(ONE/THREE) * 4**(TWO/THREE)

# Dimension = 3
RS_FACTOR = (THREE/(FOUR*pi))**(ONE/THREE)
LDA_X_FACTOR = -X_FACTOR_C

params_a_alpha = 1
lda_x_ax = -params_a_alpha*RS_FACTOR*X_FACTOR_C/2**(FOUR/THREE)

def f_lda_x(r_s, z):
    return lda_x_ax*( (1 + z)**(FOUR/THREE) + (1 - z)**(FOUR/THREE) )/r_s

def f(r_s, z):
    return f_lda_x(r_s, z)


ζ = symbols("zeta")
ρ = symbols("rho")

#r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)
r_s = symbols("r_s")

print()
print("Spin pol")
pprint(f(r_s, ζ))

# Non spinpol
ζ = 0
print()
print("Non-spin pol")
pprint(f(r_s, ζ))

ε_x = f(r_s, ζ)

print("Derivative w.r.t r_s:")
d_ε_x = diff(ε_x, r_s) 
pprint(d_ε_x)

from sympy.utilities.codegen import codegen
code1 = codegen( ("eps_x", ε_x), language="julia")
print(code1[0][1])
code1 = codegen( ("d_eps_x", d_ε_x), language="julia")
print(code1[0][1])


