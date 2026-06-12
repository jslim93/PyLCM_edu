# PyLCM Performance + Ensemble Improvement Plan

> Branch `feature/performance-ensemble` (worktree `../PyLCM_perf`), off `main` (v1.0.0).
> Constraint: never add a Claude/Co-Authored-By line. Never touch `/Users/dr.cloud/PyLCM_edu`.

**Goal:** Make single runs much faster (struct-of-arrays + numba on the condensation hot
loop and the collision shuffle/efficiency helpers) and make ensembles run efficiently in
parallel across cores — without regressing validated physics.

**Profiling basis (full run, 2000 ptcl × 400 steps):** drop_condensation 42%, collision
random.shuffle 13%, E_H80 8%, ws_drops_beard 3%. Condensation is deterministic (golden-testable);
collision is stochastic (invariant + statistical tests).

**Verification strategy:**
- **Condensation (deterministic):** golden regression — snapshot a fixed-seed condensation-only
  run from the current object code, assert the SoA/njit path reproduces it within `rtol=1e-9`.
- **Collision (stochastic):** the v1.0 conservation invariants must still pass, plus a statistical
  check that ensemble-mean diagnostics from the new path match the old within Monte-Carlo tolerance.

---

## Task 1: Golden regression harness for condensation
- Create `tests/golden/` and a script `validation/make_golden.py` that runs a fixed-seed,
  condensation-only sequence (small: 200 ptcl, 300 steps) using the CURRENT object code and
  saves the final `M`, `A`, `T`, `q` arrays to `tests/golden/condensation_golden.npz`.
- `tests/test_golden_condensation.py`: re-run the same setup, assert match within `rtol=1e-9`.
- Commit golden artifact + test. This locks current behavior BEFORE any rewrite.

## Task 2: Parallel ensemble adapter (real physics)
- Add `pylcm_run.py` with `run_single_series(seed, n_ptcl, nt, ...)` adapted from
  `ensemble_comparison.py`, returning a 1-D diagnostic time series (LWC at saved steps).
- Wire `ensemble.py.run_ensemble` to it; verify parallel across cores (joblib `n_jobs=-1`).
- `tests/test_ensemble_real.py`: run 2 small members, assert shape + seed determinism.
- Commit.

## Task 3: SoA + njit condensation (golden-guarded)
- Add a struct-of-arrays container (parallel numpy arrays M/A/Ns/kappa) + converters to/from
  the object list.
- Implement `drop_condensation_soa(...)` as an `@njit` function replicating condensation.py:8
  math exactly (scalars once, per-index loop, calls the already-jitted `radius_liquid_euler`).
- Assert it passes `tests/test_golden_condensation.py` within tolerance.
- Switch the driver path to SoA condensation; re-profile; record before/after in `validation/PROFILE.md`.
- Commit only if golden passes AND it is measurably faster; else keep object path.

## Task 4: Faster collision on SoA (invariant + statistical guarded)
- Replace `random.shuffle(list)` with `np.random.permutation` over indices.
- `@njit` `E_H80` and `ws_drops_beard`.
- v1.0 conservation invariants must pass; add a statistical ensemble-mean equivalence test
  (old vs new within Monte-Carlo tolerance). Commit only if green and faster.

## Task 5: Wrap-up
- Update `validation/PROFILE.md` with end-to-end speedup; note any deferred work.
- Run full `pytest tests/ -v`. Finish the branch (PR or merge per user).
