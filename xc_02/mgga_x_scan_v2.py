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
#def mgga_exchange(func, rs, z, xs0, xs1, u0, u1, t0, t1):
#    return lda_x_spin(rs, 1)*func(xs0, u0, t0)

def mgga_exchange(func, rs, z, xs0, xs1, u0, u1, t0, t1):
    return lda_x_spin(rs, z)*func(xs0, u0, t0) + lda_x_spin(rs, -z)*func(xs1, u1, t1)

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

z = 0.0

eps_x = f(rs, z, xt, xs0, xs1, u0, u1, t0, t1)

# CODEGEN is skipped

ρ = symbols("rho")
sigma = symbols("sigma")
tau = symbols("tau")

rs_in_rho = (3/(4*pi*ρ))**(1/3)
xs0_in_sigma = sqrt(sigma/4)/((ρ/2)**(1 + 1/DIMENSIONS))
xs1_in_sigma = sqrt(sigma/4)/((ρ/2)**(1 + 1/DIMENSIONS))
t0_in_rho = (tau/2)/((ρ/2)**(1 + 2/DIMENSIONS))
t1_in_rho = (tau/2)/((ρ/2)**(1 + 2/DIMENSIONS))

dict_mainvar = {rs: rs_in_rho, xs0: xs0_in_sigma, xs1: xs1_in_sigma, t0: t0_in_rho, t1: t1_in_rho}
eps_x = eps_x.subs(dict_mainvar)

print(eps_x)

ρ_num = 1.2
sigma_num = 0.1
tau_num = 0.1

print("ρ     = %18.10f" % ρ_num)
print("sigma = %18.10f" % sigma_num)
print("tau   = %18.10f" % tau_num)

dict_num = {ρ: ρ_num, sigma: sigma_num, tau: tau_num}
print("eps_x          = %18.10f" % (eps_x.subs(dict_num)) )

d_eps_x_drho = diff(ρ*eps_x, ρ)
print("d_eps_x_drho   = %18.10f" % (d_eps_x_drho.subs(dict_num)) )

d_eps_x_dsigma = diff(ρ*eps_x, sigma)
print("d_eps_x_dsigma = %18.10f" % (d_eps_x_dsigma.subs(dict_num)) )

d_eps_x_dtau = diff(ρ*eps_x, tau)
print("d_eps_x_dtau   = %18.10f" % (d_eps_x_dtau.subs(dict_num)) )
