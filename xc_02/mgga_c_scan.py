from sympy import *
from numeric_symbols import *

X2S        = 1/(2*(6*pi**2)**(1/3))
K_FACTOR_C = 3/10*(6*pi**2)**(2/3)

def f_zeta(z):
    return ((1 + z)**(4/3) + (1 - z)**(4/3) - 2)/(2**(4/3) - 2)

def mphi(z):
    return ((1 + z)**(2/3) + (1 - z)**(2/3))/2

def tt(rs, z, xt):
    return xt/(4*2**(1/3)*mphi(z)*sqrt(rs))

# lda_c_pw_params
params_a_pp     = [0.0, 1,  1,  1]
params_a_a      = [0.0, 0.031091, 0.015545, 0.016887]
params_a_alpha1 = [0.0, 0.21370,  0.20548,  0.11125]
params_a_beta1  = [0.0, 7.5957, 14.1189, 10.357]
params_a_beta2  = [0.0, 3.5876, 6.1977, 3.6231]
params_a_beta3  = [0.0, 1.6382, 3.3662,  0.88026]
params_a_beta4  = [0.0, 0.49294, 0.62517, 0.49671]
params_a_fz20   = 1.709921

# ifdef lda_c_pw_modified_params
params_a_a      = [0.0, 0.0310907, 0.01554535, 0.0168869]
params_a_fz20   = 1.709920934161365617563962776245


# LDA PW (* Equation (10) *)
def g_aux(k, rs):
    return params_a_beta1[k]*sqrt(rs) + params_a_beta2[k]*rs + \
        params_a_beta3[k]*rs**1.5 + params_a_beta4[k]*rs**(params_a_pp[k] + 1)

def g(k, rs):
    return -2*params_a_a[k]*(1 + params_a_alpha1[k]*rs) * log(1 +  1/(2*params_a_a[k]*g_aux(k, rs)))

# LDA PW (* Equation (8) *)
# (* Attention, the function g parametrizes -alpha *)
def f_pw(rs, zeta):
    return g(1, rs) + zeta**4*f_zeta(zeta)*(g(2, rs) - g(1, rs) + \
        g(3, rs)/params_a_fz20) - f_zeta(zeta)*g(3, rs)/params_a_fz20


#$define lda_c_pw_params
#$define lda_c_pw_modified_params
#$include "lda_c_pw.mpl"

params_a_beta  = 0.06672455060314922
params_a_gamma = (1 - log(2.0))/pi**2
params_a_BB    = 1

mgamma = params_a_gamma

def mbeta(rs, t):
    return params_a_beta

BB     = params_a_BB

def tp(rs, z, xt):
    return tt(rs, z, xt)

# (* Equation (8) *)
def A(rs, z, t):
  return mbeta(rs, t)/(mgamma*(exp(-f_pw(rs, z)/(mgamma*mphi(z)**3)) - 1))

#(* Equation (7) *)
def f1(rs, z, t):
    return t**2 + BB*A(rs, z, t)*t**4

def f2(rs, z, t):
    return mbeta(rs, t)*f1(rs, z, t)/(mgamma*(1 + A(rs, z, t)*f1(rs, z, t)))

def fH(rs, z, t):
    return mgamma*mphi(z)**3*log(1 + f2(rs, z, t))

def f_pbe(rs, z, xt, xs0, xs1):
  return f_pw(rs, z) + fH(rs, z, tp(rs, z, xt))

#$include "gga_c_regtpss.mpl"

# (* in the paper we have beta_a = 0.066725 *)
beta_a = 0.066724550603149220
beta_b = 0.1
beta_c = 0.1778

#(* we redefine beta here *)
#(* this is the Hu and Langreth expression *)
def mbeta(rs, t):
    return beta_a*(1 + beta_b*rs)/(1 + beta_c*rs)

# $include "gga_c_scan_e0.mpl"
def scan_e0_g(rs, z, t):
    return (1 + 4*A(rs, z, t)*t**2)**(-1/4)

def f2(rs, z, t):
    return mbeta(rs, t)*(1 - scan_e0_g(rs, z, t))/(mgamma*A(rs, z, t))


# include "mgga_x_scan.mpl"

scan_b1c = 0.0285764
scan_b2c = 0.0889
scan_b3c = 0.125541

def scan_eclda0(rs):
    return -scan_b1c/(1 + scan_b2c*sqrt(rs) + scan_b3c*rs)

scan_chi_infty = 0.12802585262625815

def scan_g_infty(s):
    return 1/(1 + 4*scan_chi_infty*s**2)**(1/4)

# (* in the paper it is 2.3631 *)
scan_G_cnst = 2.363

def scan_Gc(z):
    return (1 - scan_G_cnst*(2**(1/3) - 1)*f_zeta(z))*(1 - z**12)

def scan_H0(rs, s):
    return scan_b1c*log(1 + (exp(-scan_eclda0(rs)/scan_b1c) - 1)*(1 - scan_g_infty(s)))

def scan_e0(rs, z, s):
    return (scan_eclda0(rs) + scan_H0(rs, s))*scan_Gc(z)

def t_total(z, ts0, ts1):
    return (ts0*((1 + z)/2)**(5/3) + ts1*((1 - z)/2)**(5/3))

def scan_alpha(z, xt, ts0, ts1):
    return (t_total(z, ts0, ts1) - xt**2/8)/(K_FACTOR_C*t_total(z, 1, 1))

#(* set parameters of f_alpha *)
params_a_c1 = 0.64
params_a_c2 = 1.5
params_a_d  = 0.7

def scan_f_alpha(a):
    func1 = exp(-params_a_c1*a/(1 - a))
    func2 = -params_a_d*exp(params_a_c2/(1 - a))
    return Piecewise( (func1, a <= 1), (func2, True) )

def scan_f(rs, z, xt, xs0, xs1, ts0, ts1):
    return f_pbe(rs, z, xt, xs0, xs1) + scan_f_alpha(scan_alpha(z, xt, ts0, ts1))*( 
        scan_e0(rs, z, X2S*2**(1/3)*xt) - f_pbe(rs, z, xt, xs0, xs1) )

def f(rs, z, xt, xs0, xs1, us0, us1, ts0, ts1):
    return scan_f(rs, z, xt, xs0, xs1, ts0, ts1)

rs = symbols("rs", positive=True)
z, xt, xs0, xs1, us0, us1, ts0, ts1 = symbols("z xt xs0 xs1 us0 us1 ts0 ts1", real=True) 

z = 0
eps_c = f(rs, z, xt, xs0, xs1, us0, us1, ts0, ts1)

d_eps_c_d_rs = diff(f(rs, z, xt, xs0, xs1, us0, us1, ts0, ts1), rs)

from sympy.utilities.codegen import codegen

code1 = codegen( ("eps_c", eps_c), language="julia")
print(code1[0][1])

code1 = codegen( ("d_eps_c_d_rs", d_eps_c_d_rs), language="julia")
print(code1[0][1])

print("Nops = ", count_ops(eps_c))
print("Nops = ", count_ops(d_eps_c_d_rs))

ρ = 1.0
r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)
print(eps_c.subs({xt: 0.0, ts0: 0.0, ts1: 0.0, rs: r_s}))