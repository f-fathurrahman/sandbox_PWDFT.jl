from sympy import *
from numeric_symbols import *
#from symbolic_symbols import *

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

#$ifdef gga_c_pbe_params
params_a_beta  = 0.06672455060314922
params_a_gamma = (1.0 - log(2.0))/pi**2
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

# (* Equation (7) *)
def f1(rs, z, t):
    return t**2 + BB*A(rs, z, t)*t**4

def f2(rs, z, t):
    return mbeta(rs, t)*f1(rs, z, t)/(mgamma*(1 + A(rs, z, t)*f1(rs, z, t)))

def fH(rs, z, t):
    return mgamma*mphi(z)**3*log(1 + f2(rs, z, t))

def f_pbe(rs, z, xt, xs0, xs1):
    return f_pw(rs, z) + fH(rs, z, tp(rs, z, xt))

def f(rs, z, xt, xs0, xs1):
    return f_pbe(rs, z, xt, xs0, xs1)


r_s, z, xt, xs0, xs1 = symbols("r_s z xt xs0 xs1")
z = 0 # no spin pol

ε_c = f(r_s, z, xt, xs0, xs1)
from sympy.utilities.codegen import codegen

code1 = codegen( ("eps_c", ε_c), language="julia")
print(code1[0][1])

d_ε_c1 = diff(ε_c, r_s)
code1 = codegen( ("d_eps_c1", d_ε_c1), language="julia")
print(code1[0][1])

d_ε_c1 = diff(ε_c, r_s)
code1 = codegen( ("d_eps_c1", d_ε_c1), language="C")
print(code1[0][1])

print("Nops = ", count_ops(ε_c))
print("Nops = ", count_ops(d_ε_c1))

#d_ε_c2 = diff(ε_c, xt)
#code1 = codegen( ("d_eps_c2", d_ε_c2), language="julia")
#print(code1[0][1])

