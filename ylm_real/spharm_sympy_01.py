from sympy import *

init_printing(use_unicode=True)

π = pi

def gen_Nlm(l,m):
    term1 = Rational((2*l + 1) * factorial(l-m), 4*factorial(l+m))
    return sqrt(term1) / sqrt(π)

def gen_Ylm_complex(l,m,θ,ϕ):
    Nlm = gen_Nlm(l,m)
    # force the factor to be Rational
    expr1 = Rational(-1,1)**m * Nlm * assoc_legendre(l,m,cos(θ)) * exp(I*m*ϕ)
    #expr1 = Nlm * assoc_legendre(l,m,cos(θ)) * exp(I*m*ϕ) # no (-1)^m factor
    # force sqrt(1-cos(theta)**2) -> sin(θ)
    return simplify(expr1.subs( sqrt(1 - cos(θ)**2), sin(θ) ))


θ, ϕ = symbols("theta phi", real=True)

lmax = 3
Ylm = {} # initialize empty dict
for l in range(0,lmax+1):
    for m in range(-l,l+1):
        Ylm[l,m] = gen_Ylm_complex(l,m,θ,ϕ)

Ylm_real = {}
for l in range(0,lmax+1):
    for m in range(-l,l+1):
        m1m = Rational(-1,1)**m
        if m < 0:
            Ylm_real[l,m] = I/sqrt(2) * simplify( Ylm[l,m] - m1m * Ylm[l,-m] )
        elif m > 0:
            Ylm_real[l,m] = 1/sqrt(2) * simplify( Ylm[l,-m] + m1m * Ylm[l,m] )
        else: # m == 0
            Ylm_real[l,m] = Ylm[l,m]



#lmax = 1,
#Rlm = { (1,0): Ylm[1,0] }
#for m in range()
#Rlm = {
#    (-1,1): (-1)**1/sqrt(2) * (Ylm[l,m] + conj)
#}