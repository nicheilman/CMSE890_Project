import numpy as np

def walls(x):
    return np.logical_or(np.isclose(x[1], 0), np.isclose(x[1], 1))

def inflow(x):
    return np.isclose(x[0], 0)

def outflow(x):
    return np.isclose(x[0], 1)
