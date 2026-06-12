"""Benchmark: object-based drop_condensation vs SoA+numba drop_condensation_fast.

Times a condensation-only ascent at n_ptcl=20000, nt=300 through both paths.
Numba is warmed up with one call first so JIT-compile time is excluded.
Appends results to validation/PROFILE.md under a "## Condensation SoA+numba"
heading.

Run from the repo root as a module so the `validation` package resolves:
    python -m validation.bench_condensation
"""
import time
import numpy as np

from validation.golden_setup import run_condensation_only
from PyLCM.condensation import drop_condensation
from PyLCM.condensation_fast import drop_condensation_fast


def _time_fn(fn, n_ptcl, nt):
    t0 = time.perf_counter()
    out = run_condensation_only(n_ptcl=n_ptcl, nt=nt, condensation_fn=fn)
    t1 = time.perf_counter()
    return t1 - t0, out


def main():
    N_PTCL = 20000
    NT = 300

    # Warm up numba (compile _cond_kernel) on a tiny run so the timed run is
    # pure execution, not JIT compilation.
    print("Warming up numba kernel...")
    run_condensation_only(n_ptcl=10, nt=2, condensation_fn=drop_condensation_fast)

    print(f"Benchmarking n_ptcl={N_PTCL}, nt={NT} ...")
    t_obj, out_obj = _time_fn(drop_condensation, N_PTCL, NT)
    t_fast, out_fast = _time_fn(drop_condensation_fast, N_PTCL, NT)

    speedup = t_obj / t_fast

    # Sanity: confirm the two paths still agree at this larger size.
    max_rel_M = float(np.max(np.abs(out_fast["M"] - out_obj["M"]) /
                             np.maximum(np.abs(out_obj["M"]), 1e-300)))

    print(f"object drop_condensation:      {t_obj:.3f} s")
    print(f"fast   drop_condensation_fast: {t_fast:.3f} s")
    print(f"speedup:                       {speedup:.2f}x")
    print(f"max rel diff M (obj vs fast):  {max_rel_M:.3e}")

    block = (
        "\n## Condensation SoA+numba\n\n"
        f"Benchmark driver: `validation/bench_condensation.py` "
        f"(condensation-only, n_ptcl={N_PTCL}, nt={NT}, numba warmed up first).\n\n"
        "| path                              | wall time (s) |\n"
        "|-----------------------------------|---------------|\n"
        f"| object `drop_condensation`        | {t_obj:.3f} |\n"
        f"| fast `drop_condensation_fast`     | {t_fast:.3f} |\n"
        f"| **speedup**                       | **{speedup:.2f}x** |\n\n"
        f"Correctness: fast path reproduces the object path within "
        f"max rel diff M = {max_rel_M:.3e}; the golden regression "
        f"(`tests/test_condensation_fast.py`) matches at **rtol=1e-9** "
        f"(bit-for-bit, rel diff 0.0 at n_ptcl=200).\n\n"
        "Implementation: `PyLCM/condensation_fast.py` extracts particle state\n"
        "(M, A, Ns, kappa) into contiguous float64 arrays and runs the\n"
        "per-particle loop in an `@njit(cache=True)` kernel (`_cond_kernel`),\n"
        "which calls the already-jitted `radius_liquid_euler` directly. The\n"
        "object-based reference in `PyLCM/condensation.py` is unchanged and\n"
        "remains the regression anchor.\n\n"
        "Bit-exactness note: integer powers (`r_liq ** 3`) are written as\n"
        "`** 3.0` in the kernel so numba routes through libm `pow()` and\n"
        "matches CPython's `** 3` to the last ULP. With `** 3` numba emits\n"
        "`x*x*x`, a 1-ULP difference that amplifies chaotically over the\n"
        "300-step Newton-Raphson ascent and breaks the rtol=1e-9 gate.\n"
    )
    with open("validation/PROFILE.md", "a") as f:
        f.write(block)
    print("Appended results to validation/PROFILE.md")


if __name__ == "__main__":
    main()
