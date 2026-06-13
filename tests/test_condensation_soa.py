"""Persistent struct-of-arrays condensation must be bit-identical to the object
golden — proving the SoA layout is a pure performance optimization of OUR
validated physics, not a different (regressed) implementation.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")

from PyLCM.parameters import p0, r_a, cp, rv, rho_aero, z_env
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import condense_soa

GOLDEN = "tests/golden/condensation_golden.npz"


def _run_condensation_soa(n_ptcl=200, nt=300):
    """Same setup as validation/golden_setup.run_condensation_only, but with the
    particle state held as PERSISTENT numpy arrays and grown via condense_soa."""
    np.random.seed(12345)
    T0 = 293.2; P0 = 1013.0e2; RH = 0.95; w = 1.0; z0 = 0.0; zmax = 3000.0; dt = 1.0
    N_aero = np.array([118.0, 11.0, 0.72]) * 1e6
    mu_aero = np.log(np.array([0.019, 0.056, 0.46]) * 1e-6)
    sigma_aero = np.log(np.array([3.3, 1.6, 2.2]))
    kappa_arr = [1.6] * 4
    theta_profiles = T0 * (p0 / P0) ** (r_a / cp) + 5e-3 * z_env
    es = esatw(T0)
    q0 = RH * es / (P0 - RH * es) * r_a / rv

    T, q, pl = aero_init("Weighting_factor", n_ptcl, P0, z0, T0, q0,
                         N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, True)
    # extract persistent arrays ONCE (this is the whole point — no per-step rebuild)
    M = np.array([p.M for p in pl], dtype=np.float64)
    A = np.array([p.A for p in pl], dtype=np.float64)
    Ns = np.array([p.Ns for p in pl], dtype=np.float64)
    kappa = np.array([p.kappa for p in pl], dtype=np.float64)

    P, z = P0, z0
    for t in range(nt):
        z, T, P = ascend_parcel(z, T, P, w, dt, (t + 1) * dt, zmax, theta_profiles, None, "linear")
        _, _, air_mass = parcel_rho(P, T)
        T, q = condense_soa(M, A, Ns, kappa, T, q, P, dt, air_mass, rho_aero,
                            kohler_activation_radius=True, switch_kappa_koehler=True)
    return {"M": np.sort(M), "A": np.sort(A), "T": float(T), "q": float(q)}


def test_soa_condensation_matches_golden():
    ref = np.load(GOLDEN)
    out = _run_condensation_soa()
    assert np.allclose(out["M"], ref["M"], rtol=1e-9, atol=0.0)
    assert np.allclose(out["A"], ref["A"], rtol=1e-9, atol=0.0)
    assert np.isclose(out["T"], float(ref["T"]), rtol=1e-9)
    assert np.isclose(out["q"], float(ref["q"]), rtol=1e-9)
