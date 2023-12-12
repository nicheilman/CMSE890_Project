import numpy as np

def walls(x):
    #if (x[2] > 0) and (x[2] < 10): 
    return np.full(x.shape[1], True)

def inflow(x):
    return np.isclose(x[2], 10)

def outflow(x):
    return np.isclose(x[2], 0)
