"""Programmatic single-run driver for PyLCM (no notebook/widgets).

Provides `run_single_series`, a stochastic single realization returning a 1-D
diagnostic time series suitable for `ensemble.run_ensemble`. Members differ
through the RNG state, which `ensemble.run_member` seeds before each call.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

from PyLCM.parameters import *
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw, drop_condensation
from PyLCM.collision import collection
from Post_process.analysis import ts_analysis


def run_single_series(n_ptcl=2000, nt=1000, collect_every=50, diagnostic="LWC"):
    """Run one full condensation+collision ascent; return a 1-D diagnostic series.

    The RNG is NOT seeded here on purpose: `ensemble.run_member` seeds globally
    before calling, so each ensemble member is an independent realization.
    Returned series is the chosen `diagnostic` sampled every `collect_every` steps.
    """
    T0 = 293.2; P0 = 1013.0e2; RH = 0.92; w = 1.0; z0 = 0.0; zmax = 3000.0; dt = 1.0
    N_raw = np.array([118.0, 11.0, 0.72])
    mu_aero = np.log(np.array([0.019, 0.056, 0.46]) * 1e-6)
    sigma_aero = np.log(np.array([3.3, 1.6, 2.2]))
    kappa_arr = [1.6] * 4
    N_aero = N_raw * 1e6
    theta_init = T0 * (p0 / P0) ** (r_a / cp)
    theta_profiles = theta_init + 5e-3 * z_env
    e_s = esatw(T0)
    q0 = RH * e_s / (P0 - RH * e_s) * r_a / rv
    rm_spec = [1e-6, 25e-6]

    T_parcel, q_parcel, particles_list = aero_init(
        "Random", n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, False)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    series = []
    for t in range(nt):
        tt = (t + 1) * dt
        z_parcel, T_parcel, P_parcel = ascend_parcel(
            z_parcel, T_parcel, P_parcel, w, dt, tt, zmax, theta_profiles, None, "linear")
        rho_p, _, air_mass = parcel_rho(P_parcel, T_parcel)
        ct = at = et = da = 0.0
        particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = drop_condensation(
            particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass, S_lst,
            rho_aero, False, ct, at, et, da, False)
        ac = au = pr = 0.0
        particles_list, ac, au, pr = collection(
            dt, particles_list, rho_p, rho_liq, P_parcel, T_parcel, ac, au, pr,
            False, z_parcel, zmax, w)

        if (t + 1) % collect_every == 0:
            sp, qa, qc, qr, na, nc, nr, _, ra, rs = ts_analysis(
                particles_list, air_mass, rm_spec, 60, n_ptcl)
            rho2, _, _ = parcel_rho(P_parcel, T_parcel)
            value = {"LWC": (qc + qr) * rho2, "Nr": nr, "Nc": nc, "qc": qc, "qr": qr}[diagnostic]
            series.append(value)

    return np.array(series)
