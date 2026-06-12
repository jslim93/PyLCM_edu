"""Determinism guard for the Tier-1 jit of E_H80 and ws_drops_beard.

These two helpers were decorated with @jit(nopython=True, cache=True) and their
lookup tables were lifted to module level. That refactor must NOT change the
numerical result for any input. To prove it, this test embeds a faithful pure
-Python copy of the *pre-jit* implementations (the ground truth) and asserts the
current (jitted) collision.E_H80 / collision.ws_drops_beard reproduce them to
rtol <= 1e-10 over a representative grid of radii / pairs.
"""
import math
import numpy as np

from PyLCM.collision import E_H80, ws_drops_beard
from PyLCM.parameters import g


# ---------------------------------------------------------------------------
# Pre-jit reference implementations (verbatim copies of the original code path)
# ---------------------------------------------------------------------------
def _E_H80_ref(r1, r2):
    r0 = np.array([6.0, 8.0, 10.0, 15.0, 20.0, 25.0,
                   30.0, 40.0, 50.0, 60.0, 70.0, 100.0,
                   150.0, 200.0, 300.0])
    rat = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25,
                    0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                    0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
                    0.90, 0.95, 1.00])
    ecoll = np.array([
        [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
        [0.003, 0.003, 0.003, 0.004, 0.005, 0.005, 0.005, 0.010, 0.100, 0.050, 0.200, 0.500, 0.770, 0.870, 0.970],
        [0.007, 0.007, 0.007, 0.008, 0.009, 0.010, 0.010, 0.070, 0.400, 0.430, 0.580, 0.790, 0.930, 0.960, 1.000],
        [0.009, 0.009, 0.009, 0.012, 0.015, 0.010, 0.020, 0.280, 0.600, 0.640, 0.750, 0.910, 0.970, 0.980, 1.000],
        [0.014, 0.014, 0.014, 0.015, 0.016, 0.030, 0.060, 0.500, 0.700, 0.770, 0.840, 0.950, 0.970, 1.000, 1.000],
        [0.017, 0.017, 0.017, 0.020, 0.022, 0.060, 0.100, 0.620, 0.780, 0.840, 0.880, 0.950, 1.000, 1.000, 1.000],
        [0.030, 0.030, 0.024, 0.022, 0.032, 0.062, 0.200, 0.680, 0.830, 0.870, 0.900, 0.950, 1.000, 1.000, 1.000],
        [0.025, 0.025, 0.025, 0.036, 0.043, 0.130, 0.270, 0.740, 0.860, 0.890, 0.920, 1.000, 1.000, 1.000, 1.000],
        [0.027, 0.027, 0.027, 0.040, 0.052, 0.200, 0.400, 0.780, 0.880, 0.900, 0.940, 1.000, 1.000, 1.000, 1.000],
        [0.030, 0.030, 0.030, 0.047, 0.064, 0.250, 0.500, 0.800, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000],
        [0.040, 0.040, 0.033, 0.037, 0.068, 0.240, 0.550, 0.800, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000],
        [0.035, 0.035, 0.035, 0.055, 0.079, 0.290, 0.580, 0.800, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000],
        [0.037, 0.037, 0.037, 0.062, 0.082, 0.290, 0.590, 0.780, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000],
        [0.037, 0.037, 0.037, 0.060, 0.080, 0.290, 0.580, 0.770, 0.890, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000],
        [0.037, 0.037, 0.037, 0.041, 0.075, 0.250, 0.540, 0.760, 0.880, 0.920, 0.950, 1.000, 1.000, 1.000, 1.000],
        [0.037, 0.037, 0.037, 0.052, 0.067, 0.250, 0.510, 0.770, 0.880, 0.930, 0.970, 1.000, 1.000, 1.000, 1.000],
        [0.037, 0.037, 0.037, 0.047, 0.057, 0.250, 0.490, 0.770, 0.890, 0.950, 1.000, 1.000, 1.000, 1.000, 1.000],
        [0.036, 0.036, 0.036, 0.042, 0.048, 0.230, 0.470, 0.780, 0.920, 1.000, 1.020, 1.000, 1.000, 1.000, 1.000],
        [0.040, 0.040, 0.035, 0.033, 0.040, 0.112, 0.450, 0.790, 1.010, 1.030, 1.040, 1.000, 1.000, 1.000, 1.000],
        [0.033, 0.033, 0.033, 0.033, 0.033, 0.119, 0.470, 0.950, 1.300, 1.700, 2.300, 1.000, 1.000, 1.000, 1.000],
        [0.027, 0.027, 0.027, 0.027, 0.027, 0.125, 0.520, 1.400, 2.300, 3.000, 4.000, 1.000, 1.000, 1.000, 1.000],
    ]).T

    rmax = max(r1, r2)
    if rmax * 1.0E6 >= r0[14]:
        ir = 15
    else:
        for k in range(15):
            if rmax * 1e6 < r0[k]:
                ir = k
                break

    rq = min(r1 / r2, r2 / r1)
    iq = int(rq * 20)
    iq = max(iq, 1)

    if ir < 15:
        if ir >= 1:
            pp = (rmax * 1.0E6 - r0[ir - 1]) / (r0[ir] - r0[ir - 1])
            qq = (rq - rat[iq - 1]) / (rat[iq] - rat[iq - 1])
            E = (1.0 - pp) * (1.0 - qq) * ecoll[ir - 1, iq - 1] + pp * (1.0 - qq) * ecoll[ir, iq - 1] \
                + qq * (1.0 - pp) * ecoll[ir - 1, iq] + pp * qq * ecoll[ir, iq]
        else:
            qq = (rq - rat[iq - 1]) / (rat[iq] - rat[iq - 1])
            E = (1.0 - qq) * ecoll[0, iq - 1] + qq * ecoll[0, iq]
    else:
        qq = (rq - rat[iq - 1]) / (rat[iq] - rat[iq - 1])
        E = min((1.0 - qq) * ecoll[14, iq - 1] + qq * ecoll[14, iq], 1.0)

    if E < 1.0E-20:
        E = 0.0
    E = max(E, 0.0)
    return E


def _sigma_air_liq_ref(tabs):
    # Mirror of the (already jitted) condensation.sigma_air_liq used by the
    # original ws_drops_beard large-drop branch. Compared via the live one to
    # avoid drift; here we just import it.
    from PyLCM.condensation import sigma_air_liq
    return sigma_air_liq(tabs)


def _ws_drops_beard_ref(radius, rho_parcel, rho_liq, p_env, T_parcel):
    b = [-0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4, -0.327815e-5]
    c = [-0.500015e1, 0.523778e1, -0.204914e1, 0.475294, -0.542819e-1, 0.238449e-2]

    eta0 = 1.818e-5
    l0 = 6.62e-8
    p0 = 1013.25
    T0 = 293.15
    rho0 = 1.292509

    diameter = max(2.0 * radius, 0.1e-6)
    eta = rho_parcel * eta0 / rho0
    l = l0 * (eta / eta0) * (p0 / p_env) * math.sqrt(T_parcel / T0)
    Cac = 1.0 + 2.5 * l / diameter

    if diameter <= 19.0e-6:
        C1 = (rho_liq - rho_parcel) * g / (18.0 * eta)
        return C1 * Cac * diameter ** 2
    elif diameter <= 1070.0e-6:
        C2 = 4.0 * rho_parcel * (rho_liq - rho_parcel) * g / (3.0 * eta ** 2)
        NDa = C2 * diameter ** 3
        XX = math.log(NDa)
        YY = sum(b[i] * XX ** i for i in range(len(b)))
        NRe = Cac * math.exp(YY)
        return eta * NRe / (rho_parcel * diameter)
    else:
        sig = _sigma_air_liq_ref(T_parcel)
        C3 = 4.0 * (rho_liq - rho_parcel) * g / (3.0 * sig)
        Bo = C3 * diameter ** 2
        NP = sig ** 3 * rho_parcel ** 2 / (eta ** 4 * (rho_liq - rho_parcel) * g)
        XX = math.log(Bo * NP ** 0.166666666)
        YY = sum(c[i] * XX ** i for i in range(len(c)))
        NRe = NP ** 0.166666666 * math.exp(YY)
        return eta * NRe / (rho_parcel * diameter)


# Representative grid spanning Stokes / Beard intermediate / large-drop regimes.
_RADII = np.array([1e-6, 3e-6, 6e-6, 9e-6, 12e-6, 20e-6, 35e-6, 55e-6,
                   80e-6, 150e-6, 250e-6, 5e-4, 1e-3, 2e-3])


def test_E_H80_unchanged_after_jit():
    before = np.array([_E_H80_ref(r1, r2) for r1 in _RADII for r2 in _RADII])
    after = np.array([E_H80(r1, r2) for r1 in _RADII for r2 in _RADII])
    assert np.allclose(before, after, rtol=1e-10, atol=0)


def test_ws_drops_beard_unchanged_after_jit():
    p_env, T_parcel, rho_parcel, rho_liq = 90000.0, 285.0, 1.1, 1000.0
    before = np.array([_ws_drops_beard_ref(r, rho_parcel, rho_liq, p_env, T_parcel)
                       for r in _RADII])
    after = np.array([ws_drops_beard(r, rho_parcel, rho_liq, p_env, T_parcel)
                      for r in _RADII])
    assert np.allclose(before, after, rtol=1e-10, atol=0)


def test_ws_drops_beard_large_drop_regime_unchanged():
    # Force the third (Bossel/large drop) branch: diameter > 1070 um.
    p_env, T_parcel, rho_parcel, rho_liq = 80000.0, 280.0, 1.0, 1000.0
    big = np.array([600e-6, 800e-6, 1.2e-3, 2.0e-3, 3.0e-3])
    before = np.array([_ws_drops_beard_ref(r, rho_parcel, rho_liq, p_env, T_parcel)
                       for r in big])
    after = np.array([ws_drops_beard(r, rho_parcel, rho_liq, p_env, T_parcel)
                      for r in big])
    assert np.allclose(before, after, rtol=1e-10, atol=0)
