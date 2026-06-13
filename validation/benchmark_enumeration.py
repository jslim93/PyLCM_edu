"""Benchmark-only O(N^2) full-pair-enumeration collision.

NOT science, NOT shipped in the engine. Its sole purpose is the fair,
apples-to-apples speed question: with the SAME collision algorithm as the
Fortran reference (evaluate EVERY pair, no Linear-Sampling-Method subsampling),
how fast is our numba/SoA implementation?

`collide_soa` uses LSM (Shima 2009): it shuffles and evaluates only n/2 candidate
pairs per step, then upscales by n*(n-1)/(2*half). That is an O(n) ALGORITHM.
Fortran enumerates all n*(n-1)/2 pairs — an O(n^2) algorithm. Comparing the two
head-to-head conflates two different speedups:

  total speedup = (algorithm: LSM vs enumeration)  x  (implementation: numba vs Fortran)

This module isolates the second factor. `_enumerate_kernel` reuses the EXACT
per-pair physics of `_collision_kernel` (Beard terminal velocity, Hall E_H80,
Straub E_S09, p_crit acceptance, the SAM updates) but loops over every i<j pair
with NO LSM upscaling — i.e. each pair is evaluated at its true probability,
exactly as a brute-force model does.

Run:  python validation/benchmark_enumeration.py
"""
import time
import numpy as np

from PyLCM.parameters import rho_liq, pi
# The enumeration kernel now lives in the engine (collision_soa) as the
# `collision_mode="enumerate"` option; the benchmark just times it.
from PyLCM.collision_soa import _enumerate_kernel, collide_soa, seed_numba_rng


def _make_state(n, radius_um=25.0, A_each=200.0, seed=0):
    """A population of identical ~radius_um droplets, all > 10um so the full
    per-pair kernel runs (not short-circuited by the 10um skip guard)."""
    np.random.seed(seed)
    r = radius_um * 1e-6
    x = 4.0 / 3.0 * pi * rho_liq * r ** 3        # single-droplet mass
    M = np.full(n, x * A_each, dtype=np.float64)
    A = np.full(n, A_each, dtype=np.float64)
    Ns = np.full(n, x * A_each * 1e-3, dtype=np.float64)
    kappa = np.full(n, 0.5, dtype=np.float64)
    return M, A, Ns, kappa


def bench(n, dt=1.0, rho=1.0, p_env=1.0e5, T=283.0):
    seed_numba_rng(0)
    args = (dt, rho, p_env, T, False, False, False, 0.0)
    n_pairs = n * (n - 1) // 2

    # warmup JIT on a tiny state (don't time compilation)
    _enumerate_kernel(*_make_state(8), *args)
    collide_soa(*_make_state(8), dt, rho, p_env, T)

    # --- enumeration (O(n^2), all pairs) ---
    M, A, Ns, ka = _make_state(n)
    t0 = time.perf_counter()
    _enumerate_kernel(M, A, Ns, ka, *args)
    t_enum = time.perf_counter() - t0

    # --- LSM (O(n), n/2 pairs) ---
    M, A, Ns, ka = _make_state(n)
    t0 = time.perf_counter()
    collide_soa(M, A, Ns, ka, dt, rho, p_env, T)
    t_lsm = time.perf_counter() - t0

    return t_enum, n_pairs, t_lsm, n // 2


if __name__ == "__main__":
    print("Collision-step cost: O(n^2) enumeration (Fortran-style) vs O(n) LSM\n")
    print(f"{'n':>7} {'enum pairs':>13} {'enum time':>11} {'enum pairs/s':>14} "
          f"{'LSM time':>10} {'algo speedup':>13}")
    for n in (500, 1000, 2000, 4000):
        t_enum, n_pairs, t_lsm, lsm_pairs = bench(n)
        algo = t_enum / t_lsm
        print(f"{n:>7} {n_pairs:>13,} {t_enum*1e3:>9.1f}ms "
              f"{n_pairs/t_enum/1e6:>11.1f}M {t_lsm*1e3:>8.2f}ms {algo:>11.0f}x")
    print("\nAlgorithmic factor (enum pairs / LSM pairs) = (n-1)/1 ~ n.")
    print("Fortran does the O(n^2) enumeration column; our numba does it at the")
    print("'enum pairs/s' rate. Reference: Fortran 10k-particle run ~25 min "
          "(memory) => order ~1e8 pairs/s incl. condensation.")
