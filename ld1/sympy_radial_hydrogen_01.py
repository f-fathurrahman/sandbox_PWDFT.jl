from sympy import *
from sympy.physics.hydrogen import R_nl, E_nl

# Apply radial Schroedinger operator to function f with
# independent variable r
# n, l: quantum number (principal and angular momentum)
def apply_radial_operator(f, r, n, l):
    expr1 = r**2 * diff(f, r, 1)
    term1 = -1/2 * 1/r**2 * diff(expr1, r)
    Vr = -1/r + l*(l + 1)/(2*r**2)
    term2 = Vr*f
    LHS = term1 + term2
    Enl = -1/(2*n**2)
    RHS = Enl*f
    return simplify(LHS), simplify(RHS)


r = symbols("r")

R10 = R_nl(1, 0, r)

R20 = R_nl(2, 0, r)
R21 = R_nl(2, 1, r)

R30 = R_nl(3, 0, r)
R31 = R_nl(3, 1, r)
R32 = R_nl(3, 2, r)

# Check normalization
res = integrate(r**2 * R10**2, (r, 0, oo))
print("res = ", res)

# Alternative
χ = R10*r # or u(r)
res = integrate(χ**2, (r, 0, oo))
print("res = ", res)

χ = R21*r
n = 2
l = 1
Enl = -1/(2*n**2)
Vr = -1/r + l*(l + 1)/(2*r**2) - Enl
expr1 = -0.5*diff(χ, r, 2) + Vr*χ
simplify(expr1) # should be zero

LHS, RHS = apply_radial_operator(R21, r, n, l)

