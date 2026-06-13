"""Persistent struct-of-arrays parcel driver.

Ties together the SoA condensation (`condense_soa`) and collision (`collide_soa`),
both of which reuse OUR validated physics, with no per-step object<->array
conversion. This is the fast engine; every result must match the object path.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

from PyLCM.parameters import (p0, r_a, cp, rv, rho_liq, rho_aero, z_env, pi,
                              activation_radius_ts, seperation_radius_ts)
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import condense_soa
from PyLCM.collision_soa import collide_soa


def _analysis(M, A, air_mass):
    """Vectorized q/N diagnostics, matching PyLCM/Post_process classification
    (liquid radius vs activation/separation thresholds)."""
    m = A > 0
    r = np.zeros_like(M)
    r[m] = (M[m] / (A[m] * 4.0 / 3.0 * pi * rho_liq)) ** (1.0 / 3.0)
    aero = m & (r <= activation_radius_ts)
    cloud = m & (r > activation_radius_ts) & (r < seperation_radius_ts)
    rain = m & (r >= seperation_radius_ts)
    qc = np.sum(M[cloud]) / air_mass * 1e3
    qr = np.sum(M[rain]) / air_mass * 1e3
    NA = np.sum(A[aero]) / air_mass / 1e6
    NC = np.sum(A[cloud]) / air_mass / 1e6
    NR = np.sum(A[rain]) / air_mass / 1e6
    return qc, qr, NA, NC, NR


def run_soa(seed=0, n_ptcl=2000, nt=1500, dt=1.0, T0=293.2, P0=1013e2, RH=0.92,
            w=1.0, N_raw=(118., 11., .72), mu_um=(.019, .056, .46),
            sig=(3.3, 1.6, 2.2), kappa=1.6, collisions=True, switch_turb=False,
            eps=0.0, collect=None):
    """One full ascent on persistent arrays. Returns (diagnostics_by_time, (M,A))."""
    if collect is None:
        collect = (nt // 3, 2 * nt // 3, nt)
    mu = np.log(np.array(mu_um) * 1e-6)
    sg = np.log(np.array(sig))
    th = T0 * (p0 / P0) ** (r_a / cp) + 5e-3 * z_env
    q0 = RH * esatw(T0) / (P0 - RH * esatw(T0)) * r_a / rv

    np.random.seed(seed)
    T, q, pl = aero_init("Random", n_ptcl, P0, 0.0, T0, q0, np.array(N_raw) * 1e6,
                         mu, sg, rho_aero, [kappa] * (len(N_raw) + 1), False)
    # extract persistent arrays ONCE
    M = np.array([p.M for p in pl], dtype=np.float64)
    A = np.array([p.A for p in pl], dtype=np.float64)
    Ns = np.array([p.Ns for p in pl], dtype=np.float64)
    ka = np.array([p.kappa for p in pl], dtype=np.float64)

    P, z = P0, 0.0
    out = {}
    for t in range(nt):
        z, T, P = ascend_parcel(z, T, P, w, dt, (t + 1) * dt, 3000.0, th, None, "linear")
        rho_p, _, air_mass = parcel_rho(P, T)
        T, q = condense_soa(M, A, Ns, ka, T, q, P, dt, air_mass, rho_aero)
        if collisions:
            M, A, Ns, ka = collide_soa(M, A, Ns, ka, dt, rho_p, P, T,
                                       switch_turb_kernel=switch_turb, epsilon_turb=eps)[:4]
        if (t + 1) in collect:
            qc, qr, NA, NC, NR = _analysis(M, A, air_mass)
            out[t + 1] = dict(T=T - 273.15, z=z, qc=qc, qr=qr, NC=NC, NR=NR, NA=NA)
    return out, (M, A)
