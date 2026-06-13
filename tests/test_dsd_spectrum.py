import math
import numpy as np
from PyLCM.parameters import rho_liq
from PyLCM.timestep_soa import dsd_spectrum


def test_dsd_spectrum_shape_and_number():
    A = np.full(100, 1.0e6)
    r = 20e-6
    M = A * 4.0 / 3.0 * math.pi * rho_liq * r ** 3
    centers, num = dsd_spectrum(M, A, air_mass=1.0e6, n_bins=40)
    assert centers.shape == (40,) and num.shape == (40,)
    assert np.all(num >= 0)
    assert np.isclose(num.sum(), A.sum() / 1.0e6 / 1e6, rtol=1e-9)
    assert 5e-6 < centers[num.argmax()] < 60e-6
