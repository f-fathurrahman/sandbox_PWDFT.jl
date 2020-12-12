from sympy import *
from numeric_symbols import *
#from symbolic_symbols import *

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
rs = symbols("rs")

#print()
#print("Spin pol")
#pprint(f(r_s, ζ))

# Non spinpol
ζ = 0
print()
print("Non-spin pol")
eps_x = f(rs, ζ)
pprint(eps_x)

ρ = 1.1
r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)
print("eps_x = %18.10f" % (eps_x.subs({rs: r_s})))

drs_drho = -6**(1/3)/(6*pi**(1/3)*ρ**(4/3))
d_eps_x_drs = diff(eps_x, rs)
Vrho = eps_x + ρ*d_eps_x_drs*drs_drho
print("Vrho = %18.10f" % (Vrho.subs({rs: r_s})))
#print("Vxc = %18.10f" % (Vxc.subs({ρ: 1.0})))

#print("Derivative w.r.t r_s:")
#d_eps_x_drs = diff(eps_x, r_s) 
#pprint(d_eps_x_drs)

#from sympy.utilities.codegen import codegen
#code1 = codegen( ("eps_eps", eps_x), language="julia")
#print(code1[0][1])
#code1 = codegen( ("d_eps_x", d_eps_x_drs), language="julia")
#print(code1[0][1])

#ρ = 1.0
#r_s_num = (3/(4*pi*ρ))**(1/3)
#print("eps_x       = %18.10f" % eps_x.subs({r_s: r_s_num}) )

#drs_drho = -0.206783496966467*(1/ρ)**(1/3) / ρ
#print("d_eps_x_drs = %18.10f" % (d_eps_x_drs.subs({r_s: r_s_num})))
#
