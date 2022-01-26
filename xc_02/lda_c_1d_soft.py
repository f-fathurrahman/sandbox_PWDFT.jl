from sympy import *


CSC_PARAMS_SOFT_PARA = [ 0.0, # idx=0 should not be used
    #7.40, 1.120, 1.890, 0.0964,  0.0250,   2.0, 3.0, 2.431, 0.0142, 2.922
    18.40, 0.0,   7.501, 0.10185, 0.012827, 2.0, 3.0, 1.511, 0.258,  4.424
]

CSC_PARAMS_SOFT_FERRO = [ 0.0, # idx=0 should not be used
    5.24, 0.0, 1.568, 0.12856, 0.003201, 2.0, 3.0, 0.0538, 1.56e-5, 2.958
]

def f_aux(a, rs):
    return -(rs + a[5]*rs**2)*log(1 + a[8]*rs + a[9]*rs**a[10]) / \
            (2*(a[1] + a[2]*rs + a[3]*rs**a[6] + a[4]*rs**a[7]))

def f(rs, z):
    return f_aux(CSC_PARAMS_SOFT_PARA, rs) + \
           ( f_aux(CSC_PARAMS_SOFT_FERRO, rs) - f_aux(CSC_PARAMS_SOFT_PARA, rs) )*z**2


ζ = symbols("zeta")
ρ = symbols("rho")
rs = symbols("rs")

#print()
#print("Spin pol")
#pprint(f(r_s, ζ))

# Non spinpol
ζ = 0
print()
print("Non-spin pol")
eps_c = f(rs, ζ)
pprint(eps_c)

ρ = 1.1
#r_s = (3/4/pi)**(1/3)*ρ**(-1/3) # 3d
#r_s = (1/2)*ρ**(-1) # 1d
r_s = 1/(2*ρ)
print("\neps_c = %18.10f" % (eps_c.subs({rs: r_s})))

#drs_drho = -6**(1/3)/(6*pi**(1/3)*ρ**(4/3))
#d_eps_c_drs = diff(eps_c, rs)
#Vrho = eps_c + ρ*d_eps_c_drs*drs_drho
#print("Vrho = %18.10f" % (Vrho.subs({rs: r_s})))

#print("Vxc = %18.10f" % (Vxc.subs({ρ: 1.0})))

#print("Derivative w.r.t r_s:")
#d_eps_c_drs = diff(eps_c, r_s) 
#pprint(d_eps_c_drs)

#from sympy.utilities.codegen import codegen
#code1 = codegen( ("eps_eps", eps_c), language="julia")
#print(code1[0][1])
#code1 = codegen( ("d_eps_c", d_eps_c_drs), language="julia")
#print(code1[0][1])

#ρ = 1.0
#r_s_num = (3/(4*pi*ρ))**(1/3)
#print("eps_c       = %18.10f" % eps_c.subs({r_s: r_s_num}) )

#drs_drho = -0.206783496966467*(1/ρ)**(1/3) / ρ
#print("d_eps_c_drs = %18.10f" % (d_eps_c_drs.subs({r_s: r_s_num})))
#
