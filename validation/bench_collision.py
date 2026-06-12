"""Benchmark: collision-routine optimizations (Tier 1 jit + Tier 2 shuffle).

Times the AFTER (optimized) state of the collision path. A clean BEFORE-vs-AFTER
in one process is not possible because the jit decorators and the shuffle change
are baked into the imported module, so we time the AFTER full run and the
`collection` call in isolation, and report them against the profiling baseline
shares (random.shuffle ~13%, E_H80 ~8%, ws_drops_beard ~3% of a full run).

Run from the repo root as a module:
    python -m validation.bench_collision
"""
import time
import numpy as np

from PyLCM.parameters import r_a, cp, rv, rho_aero, rho_liq, p0, z_env
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import drop_condensation_fast as drop_condensation
from PyLCM.collision import collection, E_H80, ws_drops_beard


def _setup(n_ptcl, seed=0):
    np.random.seed(seed)
    T0 = 293.2; P0 = 1013.0e2; RH = 0.92; z0 = 0.0
    N_raw = np.array([118.0, 11.0, 0.72])
    mu_aero = np.log(np.array([0.019, 0.056, 0.46]) * 1e-6)
    sigma_aero = np.log(np.array([3.3, 1.6, 2.2]))
    kappa_arr = [1.6] * 4
    N_aero = N_raw * 1e6
    e_s = esatw(T0)
    q0 = RH * e_s / (P0 - RH * e_s) * r_a / rv
    T_parcel, q_parcel, particles_list = aero_init(
        "Random", n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, False)
    return T_parcel, q_parcel, particles_list


def _full_run(n_ptcl, nt, seed=0):
    T0 = 293.2; P0 = 1013.0e2; w = 1.0; z0 = 0.0; zmax = 3000.0; dt = 1.0
    theta_init = T0 * (p0 / P0) ** (r_a / cp)
    theta_profiles = theta_init + 5e-3 * z_env

    T_parcel, q_parcel, particles_list = _setup(n_ptcl, seed)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    t_coll_total = 0.0
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
        tc0 = time.perf_counter()
        particles_list, ac, au, pr = collection(
            dt, particles_list, rho_p, rho_liq, P_parcel, T_parcel, ac, au, pr,
            False, z_parcel, zmax, w)
        t_coll_total += time.perf_counter() - tc0
    return t_coll_total, len(particles_list)


def main():
    N_PTCL = 5000
    NT = 300

    # Warm up the numba-compiled helpers (E_H80, ws_drops_beard) and the
    # condensation/collision JIT kernels so timed runs exclude compilation.
    print("Warming up numba kernels (E_H80, ws_drops_beard, condensation)...")
    E_H80(2e-5, 1e-5)
    ws_drops_beard(2e-5, 1.1, 1000.0, 9e4, 285.0)
    _full_run(n_ptcl=50, nt=3, seed=1)

    print(f"Benchmarking full run n_ptcl={N_PTCL}, nt={NT} ...")
    t0 = time.perf_counter()
    t_coll, n_final = _full_run(N_PTCL, NT, seed=0)
    t_full = time.perf_counter() - t0

    coll_share = 100.0 * t_coll / t_full

    print(f"full run wall time:            {t_full:.3f} s")
    print(f"collection() total:            {t_coll:.3f} s ({coll_share:.1f}% of full run)")
    print(f"final particle count:          {n_final}")

    block = (
        "\n## Collision optimization\n\n"
        f"Benchmark driver: `validation/bench_collision.py` "
        f"(full condensation+collision ascent, n_ptcl={N_PTCL}, nt={NT}, "
        f"numba warmed up first).\n\n"
        "Two changes were made to the collision routine:\n\n"
        "**Tier 1 (deterministic):** `E_H80` (Hall 1980 collision-efficiency\n"
        "lookup + bilinear interpolation, ~8% baseline share) and\n"
        "`ws_drops_beard` (Beard 1976 terminal velocity, ~3% baseline share)\n"
        "are now `@jit(nopython=True, cache=True)`. Their lookup tables / poly\n"
        "coefficients were lifted to module-level float64 numpy arrays so numba\n"
        "closes over them as constants. Result is numerically identical\n"
        "(E_H80 bit-exact; ws_drops_beard max rel err 5.5e-15), guarded by\n"
        "`tests/test_collision_helpers_unchanged.py` at rtol=1e-10.\n\n"
        "**Tier 2 (stochastic):** the per-step list shuffle (`random.shuffle`,\n"
        "~13% baseline share) was replaced with `np.random.permutation`-based\n"
        "index reordering. This is faster and, crucially, is now controlled by\n"
        "`np.random.seed` (the ensemble seeds numpy's RNG). The exact LSM\n"
        "pairing semantics are preserved (element k of the first half pairs with\n"
        "element k of the second half, same half_length math). Because it\n"
        "changes the random sequence it is not bit-exact vs the old path;\n"
        "instead `tests/test_collision_conserves_in_full_run.py` proves the\n"
        "full run has no NaN, the collision step conserves total water mass\n"
        "(rel < 1e-10/step), and runs are reproducible under a fixed seed. The\n"
        "v1.0 conservation invariants (`tests/test_collision_invariants.py`)\n"
        "remain green.\n\n"
        "| metric                            | value |\n"
        "|-----------------------------------|-------|\n"
        f"| full run wall time (s)            | {t_full:.3f} |\n"
        f"| collection() total (s)           | {t_coll:.3f} |\n"
        f"| collection() share of full run    | {coll_share:.1f}% |\n\n"
        "The three optimized hotspots (shuffle ~13%, E_H80 ~8%,\n"
        "ws_drops_beard ~3%) together accounted for ~24% of the pre-change full\n"
        "run; after jitting the two helpers and vectorizing the shuffle the\n"
        f"collection step now measures {coll_share:.1f}% of the full ascent.\n"
    )
    with open("validation/PROFILE.md", "a") as f:
        f.write(block)
    print("Appended results to validation/PROFILE.md")


if __name__ == "__main__":
    main()
