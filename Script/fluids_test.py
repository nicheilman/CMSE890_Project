import fluids, pytest, dolfinx
import numpy as np

v = dolfinx.fem.Constant(0.0,[1,1,1])
@pytest.mark.parametrize("velocity, e", [(v, 0)])
def test_epsilon(velocity, e):
    output = fluids.epsilon(velocity)
    assert output == e

@pytest.mark.parametrize("time, flow", [(0, 0)])
def test_Inflow(time, flow):
    output = fluids.Inflow(time)
    assert output == flow

@pytest.mark.parametrize("time, pressure", [(0, -1/np.sqrt(2))])
def test_BackPressure(time, pressure):
    output = fluids.BackPressure(time)
    assert output == pressure


