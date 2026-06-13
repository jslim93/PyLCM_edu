"""Persistent struct-of-arrays parcel driver.

Ties together the SoA condensation (`condense_soa`) and collision (`collide_soa`),
both of which reuse OUR validated physics, with no per-step object<->array
conversion. This is the fast engine; every result must match the object path.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

from PyLCM.parameters import (p0, r_a, cp, rv, l_v, rho_liq, rho_aero, z_env, pi,
                              activation_radius_ts, seperation_radius_ts)
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import condense_soa
from PyLCM.collision_soa import collide_soa, seed_numba_rng


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
    qa = np.sum(M[aero]) / air_mass * 1e3
    NA = np.sum(A[aero]) / air_mass / 1e6
    NC = np.sum(A[cloud]) / air_mass / 1e6
    NR = np.sum(A[rain]) / air_mass / 1e6
    # number-weighted mean radius of cloud+rain droplets (µm), for the DSD overlay
    big = cloud | rain
    rv_mean = (np.sum(A[big] * r[big]) / np.sum(A[big]) * 1e6) if big.any() else 0.0
    return qc, qr, qa, NA, NC, NR, rv_mean


def dsd_spectrum(M, A, air_mass, n_bins=40, r_min=1e-7, r_max=5e-3):
    """Number-concentration droplet size distribution (per cm^3) over log-radius
    bins. Pure diagnostic — no physics. Returns (bin_centers_m, number_per_bin_cm3)."""
    m = A > 0
    r = np.zeros_like(M)
    r[m] = (M[m] / (A[m] * 4.0 / 3.0 * pi * rho_liq)) ** (1.0 / 3.0)
    edges = np.logspace(np.log10(r_min), np.log10(r_max), n_bins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])
    num, _ = np.histogram(r[m], bins=edges, weights=A[m])
    num = num / air_mass / 1e6
    return centers, num


def run_soa(seed=0, n_ptcl=2000, nt=1500, dt=1.0, T0=293.2, P0=1013e2, RH=0.92,
            w=1.0, N_raw=(118., 11., .72), mu_um=(.019, .056, .46),
            sig=(3.3, 1.6, 2.2), kappa=1.6, ascending_mode="linear",
            collisions=True, switch_turb=False,
            eps=0.0, lambda_ent=0.0, ihmd=0.0, init_mode="Random",
            collect=None):
    """One full ascent on persistent arrays. Returns (diagnostics_by_time, (M,A)).

    Entrainment mixing (warm-cloud, Lim & Hoffmann 2023): with lambda_ent>0 the
    parcel entrains environmental air each step and redistributes cloud liquid by
    the Inhomogeneous Mixing Degree ihmd (0 homogeneous .. 1 inhomogeneous),
    `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD` — vectorized mirror of ParameterizedMixing.
    """
    if collect is None:
        collect = (nt // 3, 2 * nt // 3, nt)
    mu = np.log(np.array(mu_um) * 1e-6)
    sg = np.log(np.array(sig))
    th = T0 * (p0 / P0) ** (r_a / cp) + 5e-3 * z_env
    q0 = RH * esatw(T0) / (P0 - RH * esatw(T0)) * r_a / rv
    # environmental vapor profile (decreases to ~2 g/kg at the top), for entrainment
    qv_prof = np.maximum(q0 - (q0 - 2e-3) / len(z_env) * np.arange(len(z_env)), 2e-3)

    # per-mode hygroscopicity: kappa may be a scalar (all modes) or a per-mode tuple
    if np.isscalar(kappa):
        k_aero = [kappa] * (len(N_raw) + 1)
    else:
        k_aero = list(kappa) + [list(kappa)[-1]]

    np.random.seed(seed)
    seed_numba_rng(seed)  # the @njit collision kernel uses Numba's separate RNG
    T, q, pl = aero_init(init_mode, n_ptcl, P0, 0.0, T0, q0, np.array(N_raw) * 1e6,
                         mu, sg, rho_aero, k_aero, False)
    # extract persistent arrays ONCE
    M = np.array([p.M for p in pl], dtype=np.float64)
    A = np.array([p.A for p in pl], dtype=np.float64)
    Ns = np.array([p.Ns for p in pl], dtype=np.float64)
    ka = np.array([p.kappa for p in pl], dtype=np.float64)

    P, z = P0, 0.0
    out = {}
    for t in range(nt):
        z, T, P = ascend_parcel(z, T, P, w, dt, (t + 1) * dt, 3000.0, th, 1200.0, ascending_mode)
        rho_p, _, air_mass = parcel_rho(P, T)
        if lambda_ent > 0.0:
            # entrainment mixing FIRST (mirror ParameterizedMixing on arrays)
            frac = min(lambda_ent * w * dt, 0.999)
            T_env = float(np.interp(z, z_env, th)) * (P / p0) ** (r_a / cp)
            q_env = float(np.interp(z, z_env, qv_prof))
            T = T + frac * (T_env - T)
            q = q + frac * (q_env - q)
            m0 = M.sum()
            M = M * (1.0 - frac)
            A = np.round(A * (1.0 - frac) ** ihmd)      # integer droplet removal
            evap = m0 - M.sum()
            q = q + evap / air_mass
            T = T - l_v * evap / cp / air_mass
        T, q = condense_soa(M, A, Ns, ka, T, q, P, dt, air_mass, rho_aero)
        if collisions:
            M, A, Ns, ka = collide_soa(M, A, Ns, ka, dt, rho_p, P, T,
                                       switch_turb_kernel=switch_turb, epsilon_turb=eps)[:4]
        if (t + 1) in collect:
            qc, qr, qa, NA, NC, NR, rv_mean = _analysis(M, A, air_mass)
            centers, num = dsd_spectrum(M, A, air_mass)
            e_s = esatw(T); e_a = q * P / (q + r_a / rv)
            out[t + 1] = dict(T=T - 273.15, T_K=T, z=z, RH=e_a / e_s, qv=q * 1e3,
                              qa=qa, qc=qc, qr=qr, NA=NA, NC=NC, NR=NR, rv=rv_mean,
                              dsd_r=centers, dsd_n=num)
    return out, (M, A)
