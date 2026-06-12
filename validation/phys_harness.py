"""Shared physics-validation harness: run a configurable parcel ascent and return
diagnostics at selected times. Used by the cross-feature physical-correctness checks.

All physics knobs are explicit so validators can vary one thing at a time.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

from PyLCM.parameters import *
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import drop_condensation_fast as drop_condensation
from PyLCM.collision import collection
from Post_process.analysis import ts_analysis


# Aerosol presets (N in cm^-3, mode radius in micron, geometric sigma)
PRESETS = {
    "maritime":    dict(N_raw=(100.0, 20.0),        mu_um=(0.08, 0.4),       sig=(1.6, 2.0),  kappa=1.0),
    "continental": dict(N_raw=(3200.0, 2900.0, 0.3), mu_um=(0.012, 0.04, 0.4), sig=(1.7, 2.0, 2.4), kappa=0.3),
    "arctic":      dict(N_raw=(15.0, 5.0),          mu_um=(0.05, 0.2),       sig=(1.6, 2.0),  kappa=0.5),
    "default":     dict(N_raw=(118.0, 11.0, 0.72),  mu_um=(0.019, 0.056, 0.46), sig=(3.3, 1.6, 2.2), kappa=1.6),
}


def run(seed=0, n_ptcl=2000, nt=1500, dt=1.0, T0=293.2, P0=1013e2, RH=0.92, w=1.0,
        aerosol="default", collisions=True, switch_turb=False, eps=0.0,
        mixing=None, collect=None):
    """Run one ascent. Returns (diagnostics_by_time, final_particles_list).

    diagnostics: dict {t: {T(C), z, qc, qr, NC, NR, NA, rv, rstd}} where NC/NR/NA
    are number concentrations (cm^-3), qc/qr are mixing ratios, rv mean volume radius.
    """
    a = PRESETS[aerosol]
    N_raw, mu_um, sig, kappa = a["N_raw"], a["mu_um"], a["sig"], a["kappa"]
    if collect is None:
        collect = (nt // 3, 2 * nt // 3, nt)

    mu = np.log(np.array(mu_um) * 1e-6)
    sg = np.log(np.array(sig))
    n_modes = len(N_raw)
    th = T0 * (p0 / P0) ** (r_a / cp) + 5e-3 * z_env
    q0 = RH * esatw(T0) / (P0 - RH * esatw(T0)) * r_a / rv

    np.random.seed(seed)
    T, q, pl = aero_init("Random", n_ptcl, P0, 0.0, T0, q0,
                         np.array(N_raw) * 1e6, mu, sg, rho_aero,
                         [kappa] * (n_modes + 1), False)
    P = P0; z = 0.0; S = 0.0
    out = {}
    for t in range(nt):
        z, T, P = ascend_parcel(z, T, P, w, dt, (t + 1) * dt, 3000.0, th, None, "linear")
        rp, _, am = parcel_rho(P, T)
        if mixing is not None:
            pl, T, q = mixing.apply(pl, T, q, P, z, dt, w, am)
        c = a_ = e = d = 0.0
        pl, T, q, S, c, a_, e, d = drop_condensation(
            pl, T, q, P, nt, dt, am, S, rho_aero, False, c, a_, e, d, False)
        if collisions:
            ac = au = pr = 0.0
            pl, ac, au, pr = collection(
                dt, pl, rp, rho_liq, P, T, ac, au, pr, False, z, 3000.0, w,
                False, False, switch_turb, eps)
        if (t + 1) in collect:
            spec, qa, qc, qr, NA, NC, NR, _, rvol, rstd = ts_analysis(
                pl, am, [1e-6, 25e-6], 60, n_ptcl)
            out[t + 1] = dict(T=T - 273.15, z=z, qc=qc, qr=qr,
                              NC=NC, NR=NR, NA=NA, rv=rvol, rstd=rstd)
    return out, pl


def total_water(pl, q, air_mass):
    """Total water mass (kg): vapor + liquid."""
    return q * air_mass + sum(p.M for p in pl)
