"""Self-contained profiling driver for PyLCM condensation hot path.

Runs ONE short condensation-only simulation directly (no widgets, no full
ensemble) so it can be wrapped in cProfile to find the dominant pure-numeric
helpers. Keep total runtime well under ~30 s.

Call path mirrors PyLCM/timestep_routine.py:
    aero_init(...) -> (T, q, particles_list)
    loop: drop_condensation(...) with P held constant by the caller.

Usage:
    python -m cProfile -s cumtime validation/profile_run.py 2>/dev/null | head -30
    python validation/profile_run.py            # plain wall-time run
"""
import time
import numpy as np

from PyLCM.aero_init import aero_init
from PyLCM.parcel import parcel_rho
from PyLCM.condensation import drop_condensation

# ----- configuration (kept small so the run finishes quickly) -----
N_PTCL = 200
N_STEPS = 3000
DT = 0.5  # s

# Parcel initial state (held constant during the run; this is a pure
# condensation micro-physics stress test, not a full ascending parcel).
P_PARCEL = 100000.0   # Pa
Z_PARCEL = 0.0        # m
T_PARCEL = 290.0      # K

# Aerosol: single accumulation mode. Convention (see MEMORY / aero_init):
#   mu = log(microns * 1e-6) ,  sigma = log(geometric_std)
N_aero = np.array([100.0e6])                 # m^-3
mu_aero = np.array([np.log(0.05 * 1e-6)])    # 0.05 micron median radius
sigma_aero = np.array([np.log(1.6)])         # geometric std 1.6
rho_aero = 1777.0                            # kg/m^3 (ammonium sulfate-ish)
k_aero = np.array([0.61])                    # hygroscopicity (kappa)
switch_kappa_koehler = True

# RH -> q : start slightly subsaturated so growth/evaporation both exercise.
RH = 0.98
from PyLCM.condensation import esatw
from PyLCM.parameters import r_a, rv
e_s = esatw(T_PARCEL)
q_parcel = RH * e_s / (P_PARCEL - RH * e_s) * r_a / rv


def build():
    T, q, particles_list = aero_init(
        "Weighting_factor", N_PTCL, P_PARCEL, Z_PARCEL, T_PARCEL, q_parcel,
        N_aero, mu_aero, sigma_aero, rho_aero, k_aero, switch_kappa_koehler,
    )
    return T, q, particles_list


def run(T, q, particles_list, n_steps=N_STEPS):
    _, _, air_mass_parcel = parcel_rho(P_PARCEL, T)
    S_lst = 0.0
    con_ts = act_ts = evp_ts = dea_ts = 0.0
    kohler_activation_radius = True
    for _ in range(n_steps):
        (particles_list, T, q, S_lst,
         con_ts, act_ts, evp_ts, dea_ts) = drop_condensation(
            particles_list, T, q, P_PARCEL, n_steps, DT, air_mass_parcel,
            S_lst, rho_aero, kohler_activation_radius,
            con_ts, act_ts, evp_ts, dea_ts, switch_kappa_koehler,
        )
    return T, q


if __name__ == "__main__":
    T, q, particles_list = build()
    t0 = time.perf_counter()
    T_out, q_out = run(T, q, particles_list)
    elapsed = time.perf_counter() - t0
    print(f"n_ptcl={N_PTCL} steps={N_STEPS} dt={DT}")
    print(f"final T={T_out:.4f} K  q={q_out:.6e}")
    print(f"wall time: {elapsed:.3f} s")
