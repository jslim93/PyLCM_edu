import numpy as np
from PyLCM.aero_init import aero_init


def _init(mode):
    mu = np.log(np.array([0.02e-6, 0.2e-6])); sigma = np.log(np.array([1.4, 1.6]))
    return aero_init(mode, 200, 95000.0, 0.0, 285.0, 0.008,
                     np.array([100e6, 20e6]), mu, sigma, 1777.0,
                     np.array([0.5, 0.5]), True)[2]


def test_random_init_multiplicities_are_integers_ge_1():
    pl = _init("Random")
    assert pl and all(p.A == round(p.A) and p.A >= 1 for p in pl)


def test_weighting_factor_init_multiplicities_are_integers_ge_1():
    pl = _init("Weighting_factor")
    assert pl and all(p.A == round(p.A) and p.A >= 1 for p in pl)
