import numpy as np
import matplotlib; matplotlib.use("Agg")
from PyLCM import parameters
from PyLCM.parcel import parcel_rho


def test_air_mass_constant_drives_parcel_rho():
    rho, V, air_mass = parcel_rho(95000.0, 285.0)
    assert abs(air_mass - parameters.PARCEL_AIR_MASS) / parameters.PARCEL_AIR_MASS < 1e-9
    assert abs(V - parameters.PARCEL_AIR_MASS / rho) / V < 1e-9
