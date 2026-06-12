import numpy as np
from PyLCM.aero_init import r_equi


def test_equilibrium_radius_positive():
    r = r_equi(-0.01, 280.0, 5.0e-8, 1777.0, True, 0.5)
    assert np.isfinite(r) and r > 0


def test_equilibrium_radius_grows_with_humidity():
    r_dry = r_equi(-0.05, 280.0, 5.0e-8, 1777.0, True, 0.5)
    r_humid = r_equi(-0.005, 280.0, 5.0e-8, 1777.0, True, 0.5)
    assert r_humid >= r_dry
