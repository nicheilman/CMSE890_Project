from ufl import (Identity, div, dot, ds, dx, inner, lhs, nabla_grad, rhs, sym)
import numpy as np

# Define strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p, den):
    return 2 * den * epsilon(u) - p * Identity(len(u))

# Define Inflow Time-Dependance
def Inflow(t):
    return (np.sin(2*np.pi*t) - np.cos(4*np.pi*t) + 1)

def BackPressure(x):
    return (np.sin(2*np.pi*x)**2 + np.sin(2*np.pi*(x - 0.125)))
