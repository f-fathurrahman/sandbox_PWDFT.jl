from sympy import *
from numeric_symbols import *
#from symbolic_symbols import *

X2S        = 1/(2*(6*pi**2)**(1/3))
K_FACTOR_C = 3/10*(6*pi**2)**(2/3)
MU_GE      = 10/81

# from the C source
params_a_c1 = 0.667
params_a_c2 = 0.8
params_a_d = 1.24
params_a_k1 = 0.065

X_FACTOR_C = THREE/EIGHT*(THREE/pi)**(ONE/THREE) * 4**(TWO/THREE)
# Dimension = 3
RS_FACTOR = (THREE/(FOUR*pi))**(ONE/THREE)
LDA_X_FACTOR = -X_FACTOR_C
DIMENSIONS = 3

def lda_x_spin(rs, z):
    return LDA_X_FACTOR*((1 + z)/2)**(1 + 1/DIMENSIONS)*(RS_FACTOR/rs)

# polarization: ferr
def mgga_exchange(func, rs, z, xs0, xs1, u0, u1, t0, t1):
    return lda_x_spin(rs, 1)*func(xs0, u0, t0)

#def mgga_exchange(func, rs, z, xs0, xs1, u0, u1, t0, t1):
#    return lda_x_spin(rs, z)*func(xs0, u0, t0) + lda_x_spin(rs, -z)*func(xs1, u1, t1)

def scan_p(x):
    return X2S**2 * x**2

def scan_alpha(x, t):
    return (t - x**2/8)/K_FACTOR_C

def scan_f_alpha(a):
    func1 = exp(-params_a_c1*a/(1 - a))
    func2 = -params_a_d*exp(params_a_c2/(1 - a))
    return Piecewise( (func1, a <= 1), (func2, True) )

def scan_h1x(x):
    return 1 + params_a_k1*(1 - params_a_k1/(params_a_k1 + x))

scan_b2 = sqrt(5913/405000)
scan_b1 = (511/13500)/(2*scan_b2)
scan_b3 = 1/2
scan_b4 = MU_GE**2/params_a_k1 - 1606/18225 - scan_b1**2

def scan_y(x, a):
    return MU_GE*scan_p(x) + scan_b4*scan_p(x)**2*exp(-scan_b4*scan_p(x)/MU_GE) + \
      (scan_b1*scan_p(x) + scan_b2*(1 - a)*exp(-scan_b3*(1 - a)**2))**2

scan_a1 = 4.9479
def scan_gx(x):
    return 1 - exp(-scan_a1/sqrt(X2S*x))

scan_h0x = 1.174
# u is not used
def scan_f(x, u, t):
    return (scan_h1x(scan_y(x, scan_alpha(x, t)))*(1 - scan_f_alpha(scan_alpha(x, t))) + 
        scan_h0x*scan_f_alpha(scan_alpha(x, t)))*scan_gx(x)

# xt is not used
def f(rs, z, xt, xs0, xs1, u0, u1, t0, t1):
    return mgga_exchange(scan_f, rs, z, xs0, xs1, u0, u1, t0, t1)

rs = symbols("rs", positive=True)
z, xt, xs0, xs1, u0, u1, t0, t1 = symbols("z xt xs0 xs1 u0 u1 t0 t1", real=True) 

z = 0

eps_x = f(rs, z, xt, xs0, xs1, u0, u1, t0, t1)
d_eps_x_d_rs = diff(eps_x, rs)
d_eps_x_d_t0 = diff(eps_x, t0)
d_eps_x_d_xs0 = diff(eps_x, xs0)

from sympy.utilities.codegen import codegen

#code1 = codegen( ("eps_x", eps_x), language="julia")
#print(code1[0][1])
#
#code1 = codegen( ("d_eps_x_d_rs", d_eps_x_d_rs), language="julia")
#print(code1[0][1])
#
#code1 = codegen( ("d_eps_x_d_t0", d_eps_x_d_t0), language="julia")
#print(code1[0][1])
#
#code1 = codegen( ("d_eps_x_d_xs0", d_eps_x_d_xs0), language="julia")
#print(code1[0][1])

print("Nops = ", count_ops(eps_x))
#print("Nops = ", count_ops(d_eps_x_d_rs))
#print("Nops = ", count_ops(d_eps_x_d_t0))
#print("Nops = ", count_ops(d_eps_x_d_xs0))


# rs and ρ relationship
ρ = symbols("rho")
r_s = (3/(4*pi*ρ))**(1/3)
#pprint(r_s)
#print("\nr_s and ρ:")
#pprint( diff(r_s, ρ) )

# s and ∇ρ relationship
delρ = symbols("∇ρ")
#s = abs(delρ)/( 2*(3*pi**2)**(1/3) * ρ**4/3 )
s = delρ/( 2*(3*pi**2)**(1/3) * ρ**4/3 )
#pprint(s)
#print("\ns and delρ:")
#pprint( diff(s, delρ) )

#rs_num = r_s.subs({ρ: 1.0})
#print("rs_num = ", rs_num)
#print("eps_x  = ", eps_x.subs({rs: rs_num, xs0: 0.0, t0: 0.0}))

ρ = 1.1
sigma = 0.0
s = sqrt(sigma)/(2 * (3*pi**2)**(1/3) * ρ**(4/3) )

r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)

print("eps_x = %18.10f" % (eps_x.subs({xs0: 0.0, xs1: 0.0, t0: 0.0, t1: 0.0, rs: r_s})) )