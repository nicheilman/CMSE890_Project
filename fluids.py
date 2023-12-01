from ufl import (FacetNormal, FiniteElement, Identity, TestFunction, TrialFunction, VectorElement,div, dot, ds, dx, inner, lhs, nabla_grad, rhs, sym)

# Define strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p, den):
    return 2 * den * epsilon(u) - p * Identity(len(u))

