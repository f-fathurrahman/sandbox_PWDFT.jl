from sympy import *
from symbolic_symbols import *

ρ = symbols("rho")
sigma = symbols("sigma")
tau = symbols("tau")

DIMENSIONS = THREE

r_s = (THREE/FOUR/pi)**(ONE/THREE)*ρ**(-ONE/THREE)
print(diff(r_s, ρ))

xs0 = sqrt(sigma/4)/((ρ/2)**(1 + 1/DIMENSIONS))
xs1 = sqrt(sigma/4)/((ρ/2)**(1 + 1/DIMENSIONS))
print(diff(xs0, sigma))
print(diff(xs1, sigma))

t0 = (tau/2)/((ρ/2)**(1 + 2/DIMENSIONS))
t1 = (tau/2)/((ρ/2)**(1 + 2/DIMENSIONS))
print(diff(t0, tau))
print(diff(t1, tau))
