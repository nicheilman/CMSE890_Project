import numpy as np

def walls(x):
    side1 = np.logical_or(np.isclose(x[1], 0), np.isclose(x[1], 1))
    side2 = np.logical_or(np.isclose(x[2], 0), np.isclose(x[2], 1))
    return np.logical_or(side1,side2)

def inflow(x):
    return np.isclose(x[0], 0)

def outflow(x):
    return np.isclose(x[0], 1)
