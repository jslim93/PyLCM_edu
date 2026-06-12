# PyLCM Condensation Profiling

Driver: `validation/profile_run.py`
Profiled with: `python -m cProfile -s cumtime validation/profile_run.py`

## Configuration

| Param            | Value                                  |
|------------------|----------------------------------------|
| n_ptcl           | 200                                    |
| steps            | 3000                                   |
| dt               | 0.5 s                                  |
| routine          | `drop_condensation` only (P constant)  |
| aerosol          | 1 accumulation mode, kappa-Köhler on   |
| machine          | Mac mini M4                            |

The driver mirrors the real call path in `PyLCM/timestep_routine.py`:
`aero_init("Weighting_factor", ...)` -> loop calling `drop_condensation`.
It is condensation-only (no ascending parcel, no collision) because
condensation is the dominant per-timestep cost.

## Baseline (BEFORE any change)

Total wall time of the `run()` loop: **~0.83-0.89 s** for 3000 steps x 200
particles (~600k particle-steps). cProfile-instrumented loop time ~0.89 s.

Top cumulative-time entries that are actual *runtime* work (import overhead
omitted):

| ncalls | tottime | cumtime | function                                  |
|--------|---------|---------|-------------------------------------------|
| 1      | 0.001   | 0.890   | `profile_run.py:run`                      |
| 3000   | 0.888   | 0.889   | `condensation.py:drop_condensation`       |

Everything else in the cProfile output (`inspect`, `sre_parse`,
`matplotlib.artist`, `numba ... dispatcher/caching`, `importlib`) is
one-time **import / JIT-compile** overhead incurred at startup, not
per-timestep cost.

### Key observation

`drop_condensation` has `tottime ≈ cumtime` (0.888 ≈ 0.889). That means
essentially all of its cost is **inlined arithmetic in its own per-particle
Python loop** — not in callees. The numeric helpers it calls
(`esatw`, `sigma_air_liq`, `radius_liquid_euler`) are **already**
`@jit(nopython=True, cache=True)` and do not appear as hot Python.

## JIT decision: NO change applied

There is **no safe pure-numeric leaf helper to JIT**:

- The remaining cost lives inside `drop_condensation`'s loop, which iterates
  over `particles` class instances and reads/writes attributes
  (`particle.M`, `.A`, `.Ns`, `.kappa`). Per the task constraints, functions
  that touch `particles` instances or Python lists of them MUST NOT be
  JIT-decorated.
- The only standalone pure-numeric helpers in the package
  (`lognormal_pdf`, `r_equi` in `aero_init.py`) run **only during
  initialization**, not in the timestep loop, so JIT-ing them yields no
  measurable runtime win and was not done.
- The already-jitted trio (`esatw`, `sigma_air_liq`, `radius_liquid_euler`)
  is the correct and only safe set; they are already optimized.

A real further speedup would require restructuring `drop_condensation` to
operate on contiguous NumPy arrays of particle state instead of a Python
list of objects (SoA refactor) so the loop body could be `nopython`-jitted.
That is a non-trivial rewrite and is explicitly **out of scope** for this
time-boxed pass. No physics was changed.

## After

N/A — no code change applied (no safe target). All 11 tests still pass
(`python -m pytest tests/ -v`); baseline correctness is unaffected because
the profiling driver is additive and does not modify any physics module.

## Condensation SoA+numba

Benchmark driver: `validation/bench_condensation.py` (condensation-only, n_ptcl=20000, nt=300, numba warmed up first).

| path                              | wall time (s) |
|-----------------------------------|---------------|
| object `drop_condensation`        | 7.425 |
| fast `drop_condensation_fast`     | 1.370 |
| **speedup**                       | **5.42x** |

Correctness: fast path reproduces the object path within max rel diff M = 0.000e+00; the golden regression (`tests/test_condensation_fast.py`) matches at **rtol=1e-9** (bit-for-bit, rel diff 0.0 at n_ptcl=200).

Implementation: `PyLCM/condensation_fast.py` extracts particle state
(M, A, Ns, kappa) into contiguous float64 arrays and runs the
per-particle loop in an `@njit(cache=True)` kernel (`_cond_kernel`),
which calls the already-jitted `radius_liquid_euler` directly. The
object-based reference in `PyLCM/condensation.py` is unchanged and
remains the regression anchor.

Bit-exactness note: integer powers (`r_liq ** 3`) are written as
`** 3.0` in the kernel so numba routes through libm `pow()` and
matches CPython's `** 3` to the last ULP. With `** 3` numba emits
`x*x*x`, a 1-ULP difference that amplifies chaotically over the
300-step Newton-Raphson ascent and breaks the rtol=1e-9 gate.
