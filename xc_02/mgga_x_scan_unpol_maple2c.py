from sympy import *

def POW_1_3(x):
    return x**(Integer(1)/Integer(3))

rho = symbols("rho")
sigma = symbols("sigma")
tau = symbols("tau")

params_a_c1 = 0.667
params_a_c2 = 0.8
params_a_d = 1.24
params_a_k1 = 0.065

M_CBRT2 = 2**(Integer(1)/Integer(3))
M_CBRT3 = 3**(Integer(1)/Integer(3))
M_CBRT4 = 4**(Integer(1)/Integer(3))
M_CBRT6 = 6**(Integer(1)/Integer(3))

t2 = M_CBRT3
t4 = POW_1_3(0.1e1 / pi)
t5 = t2 * t4
t6 = M_CBRT4
t7 = t6 * t6
t8 = t5 * t7
t9 = M_CBRT2
t10 = t9 * t9
t11 = POW_1_3(rho)
t12 = t10 * t11
t13 = M_CBRT6
t14 = pi * pi
t15 = POW_1_3(t14)
t16 = t15 * t15
t17 = 0.1e1 / t16
t18 = t13 * t17
t19 = sigma * t10
t20 = rho * rho
t21 = t11 * t11
t22 = t21 * t20
t23 = 0.1e1 / t22
t24 = t19 * t23
t25 = t18 * t24
t29 = 0.100e3 / 0.6561e4 / params_a_k1 - 0.73e2 / 0.648e3
t30 = t13 * t13
t32 = t15 * t14
t33 = 0.1e1 / t32
t34 = t29 * t30 * t33
t35 = sigma * sigma
t36 = t35 * t9
t37 = t20 * t20
t38 = t37 * rho
t40 = 0.1e1 / t11 / t38
t45 = exp(-0.27e2 / 0.80e2 * t29 * t13 * t17 * t24)
t46 = t40 * t45
t50 = sqrt(0.146e3)
t51 = t50 * t13
t52 = t51 * t17
t55 = tau * t10
t56 = t21 * rho
t57 = 0.1e1 / t56
t60 = t55 * t57 - t24 / 0.8e1
t63 = 0.5e1 / 0.9e1 * t60 * t13 * t17
t64 = 0.1e1 - t63
t66 = t64 * t64
t68 = exp(-t66 / 0.2e1)
t71 = 0.7e1 / 0.12960e5 * t52 * t24 + t50 * t64 * t68 / 0.100e3
t72 = t71 * t71
t73 = params_a_k1 + 0.5e1 / 0.972e3 * t25 + t34 * t36 * t46 / 0.288e3 + t72
t78 = 0.1e1 + params_a_k1 * (0.1e1 - params_a_k1 / t73)

t79 = t63 <= 0.1e1

t80 = params_a_c1 * t60
t81 = 0.1e1 / t64
t82 = t18 * t81
t85 = exp(-0.5e1 / 0.9e1 * t80 * t82)
t87 = exp(params_a_c2 * t81)

t89 = Piecewise( (t85, t79), (-params_a_d * t87, True) ) 

t90 = 0.1e1 - t89
t93 = t78 * t90 + 0.1174e1 * t89
t94 = sqrt(0.3e1)
t95 = 0.1e1 / t15
t96 = t30 * t95
t97 = sqrt(sigma)
t99 = t11 * rho
t100 = 0.1e1 / t99
t102 = t96 * t97 * t9 * t100
t103 = sqrt(t102)
t107 = exp(-0.98958000000000000000e1 * t94 / t103)
t108 = 0.1e1 - t107
t109 = t93 * t108
t111 = t8 * t12 * t109

eps_x = -0.3e1 / 0.16e2 * t111

print(eps_x)
print("Nops = ", count_ops(eps_x))
print(N(eps_x.subs({rho: 1.1, sigma: 0.0, tau: 0.0})))
