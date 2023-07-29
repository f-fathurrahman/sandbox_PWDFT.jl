from sympy import *
from sympy.physics.hydrogen import R_nl, E_nl

init_printing()

x = symbols("x")
r = symbols("r")

dict_subs = {
    r : exp(x)
}

R10 = R_nl(1, 0, r)

R20 = R_nl(2, 0, r)
R21 = R_nl(2, 1, r)

R30 = R_nl(3, 0, r)
R31 = R_nl(3, 1, r)
R32 = R_nl(3, 2, r)

χ = R21*r
# Function y(x)
y = (1/sqrt(r) * χ).subs(dict_subs)
n = 2
l = 1
Enl = -1/(2*n**2)
Vr = -1/r
op_V = ( 2*r**2 * (Enl - Vr) - (l + 1/2)**2 ).subs(dict_subs)
expr1 = diff(y, x, 2) + op_V*y
res = simplify(expr1) # should be zero
print("res = ", res)
