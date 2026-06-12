"""Deterministic condensation-only run, shared by the golden generator and test.

Condensation has NO random component, so with a fixed initial particle population
and a deterministic parcel ascent the final state is reproducible bit-for-bit
(within float tolerance). This is the regression anchor that guards the SoA/numba
rewrite: any fast path must reproduce this within rtol.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

from PyLCM.parameters import *
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw, drop_condensation


def run_condensation_only(n_ptcl=200, nt=300, condensation_fn=drop_condensation):
    """Run a fixed-seed, condensation-only ascent and return final state.

    `condensation_fn` lets a fast (SoA/njit) implementation be swapped in for the
    object-based reference while reusing the exact same setup and driver loop.
    Returns dict with sorted droplet masses, weights, and final T/q.
    """
    np.random.seed(12345)  # determinism even if the init mode samples

    T0 = 293.2; P0 = 1013.0e2; RH = 0.95; w = 1.0; z0 = 0.0; zmax = 3000.0; dt = 1.0
    N_raw = np.array([118.0, 11.0, 0.72])
    mu_aero = np.log(np.array([0.019, 0.056, 0.46]) * 1e-6)
    sigma_aero = np.log(np.array([3.3, 1.6, 2.2]))
    kappa_arr = [1.6] * 4
    N_aero = N_raw * 1e6
    theta_init = T0 * (p0 / P0) ** (r_a / cp)
    theta_profiles = theta_init + 5e-3 * z_env
    e_s = esatw(T0)
    q0 = RH * e_s / (P0 - RH * e_s) * r_a / rv

    T_parcel, q_parcel, particles_list = aero_init(
        "Weighting_factor", n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, True)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    for t in range(nt):
        tt = (t + 1) * dt
        z_parcel, T_parcel, P_parcel = ascend_parcel(
            z_parcel, T_parcel, P_parcel, w, dt, tt, zmax, theta_profiles, None, "linear")
        _, _, air_mass = parcel_rho(P_parcel, T_parcel)
        ct = at = et = da = 0.0
        particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = condensation_fn(
            particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass, S_lst,
            rho_aero, True, ct, at, et, da, True)

    M = np.sort(np.array([p.M for p in particles_list]))
    A = np.array([p.A for p in particles_list])
    return {"M": M, "A": np.sort(A), "T": float(T_parcel), "q": float(q_parcel)}
