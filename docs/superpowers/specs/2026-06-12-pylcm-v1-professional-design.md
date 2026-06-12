# PyLCM v1.0 "Professional" — Design Spec

**Date:** 2026-06-12
**Author:** J. Lim (jslim93) with agent assistance
**Base branch:** `feature/lecture-readiness-improvements`
**Release branch:** `release/v1.0` (in git worktree `../PyLCM_v1.0`)
**Status:** Approved design — pending implementation plan

---

## 1. Motivation

PyLCM is an educational Lagrangian Cloud Model parcel simulator. It works, but it
is not yet something an external university instructor can clone and rely on:

- **Repo hygiene** — `__pycache__/`, `.DS_Store`, `*.pyc`, `*.nbc`/`*.nbi` are tracked
  or littering the tree; no `.gitignore` discipline; no `LICENSE`.
- **Reproducibility** — no pinned environment; README tells students to `conda install`
  ~7 packages one at a time. No guarantee two students get the same stack.
- **Trust** — no tests, no CI. An external user explicitly questioned whether the
  **collision results are correct** versus the original Fortran / SAM6-LCM reference.
- **Speed** — single runs and especially ensembles are slow (≈36 min for a 100k-particle
  ensemble run); ensemble capability exists only as an ad-hoc script.
- **Placeholders** — `entrainment.py` is flagged "under development, don't use" and is
  not physically validated; students can hit broken physics silently.

The goal of v1.0 is **trustworthiness and professionalism** as an engineering release.
New unvalidated physics (a real entrainment scheme) is explicitly deferred to v1.1 so it
can be designed and validated on its own, rather than undermining the trust this release
is meant to establish.

## 2. Goals / Non-goals

**Goals**
- One-command, reproducible install (pinned env).
- Clean repo: `.gitignore`, no tracked build artifacts, a `LICENSE`.
- A pytest suite of **physics invariants** + GitHub Actions CI on every push/PR.
- A runnable **collision-validation** notebook that reproduces the PyLCM-vs-Fortran/SAM
  comparison and states the correctness story plainly.
- A proper **ensemble** module that runs members in parallel across cores, plus profiled
  single-run speedups, with measured before/after numbers.
- Entrainment (and any other non-validated stub) **hard-gated** behind an explicit
  `experimental` flag so it cannot be used silently.
- A rewritten README and a tagged `v1.0.0` release; `main` becomes the clean baseline.

**Non-goals (deferred to v1.1, separate brainstorm → spec → validation)**
- A physically-correct entrainment implementation (entrainment-rate formulation,
  homogeneous vs inhomogeneous mixing, what dilutes).
- Aggressive vectorization / algorithmic rewrites of validated condensation/collision physics.

## 3. The collision-correctness story (background that motivates the validation lane)

Source comparison across PyLCM `collision.py`, SAM6-LCM `micro_coll.f90`, and the
Fortran box model `func_coll.f90` (see `COLLISION_COMPARISON.md` and the POLARIS vault):

- **[FACT]** All three share the Linear Sampling Method (Shima et al. 2009), the Hall (1980)
  collision-efficiency table, Beard (1976) terminal velocity, and the multi-collision
  `p_crit` limiter. The core algorithm is the same.
- **[FACT]** SAM6-LCM's turbulent (Wang/Ayala) branch **omits** the Straub (2009) coalescence
  efficiency `E_S09` (`micro_coll.f90:710`); PyLCM **includes** it (`collision.py:205`).
  Where turbulence is active SAM over-counts collisions by ≈`1/E_S09` (≈1.5–10× depending
  on Weber number). PyLCM is the more physically complete code here.
- **[FACT]** SAM also omits `E_S09` in the gravitational case; PyLCM and the Fortran *box*
  model both include it — so PyLCM agrees with the box model, not with SAM.
- **[INFERRED]** Residual peak-`Nr` gap vs the box model (~3×) is LSM vs full O(N²) enumeration
  plus Monte Carlo variance — a sampling difference, not a physics error.

**Conclusion to encode in the validation notebook:** PyLCM's collision differs from SAM6
because PyLCM applies coalescence efficiency in the turbulent kernel and SAM does not;
PyLCM matches the box-model reference on the gravitational case. This is a teaching asset,
not a bug.

## 4. Worktree & branch strategy

1. On `feature/lecture-readiness-improvements`, commit **only** the good uncommitted work
   (the `PyLCM_Part1_Foundations.ipynb` / `PyLCM_Part2_Experiments.ipynb` notebooks and the
   validation/ensemble scripts) — **not** `__pycache__`/`.DS_Store`/`.pyc`/`.nbc`/`.nbi`.
2. `git worktree add ../PyLCM_v1.0 -b release/v1.0 feature/lecture-readiness-improvements`.
   The original `/Users/dr.cloud/PyLCM_edu` working dir is left untouched so the user's
   **currently-running** PyLCM session is not disturbed. All lane work happens in the worktree.
3. When tests are green, merge `release/v1.0` → `main` and tag `v1.0.0`; `main` becomes the
   professional baseline universities clone.

**Constraint:** No "Claude" coauthor / Co-Authored-By line on any commit.

## 5. Work lanes (agent-driven, parallelized by file-independence)

| Lane | Deliverables | Touches |
|------|-------------|---------|
| **A. Hygiene** | Comprehensive `.gitignore`; `git rm --cached` tracked pycache/pyc; purge `.DS_Store`/`.nbc`/`.nbi`; add `LICENSE` (MIT) | repo root |
| **B. Packaging** | Pinned `environment.yml` + `requirements.txt`; `pyproject.toml` so `pip install -e .` makes `import PyLCM` work | repo root |
| **C. Tests + validation** | `tests/` pytest (see §6); collision-validation notebook in `validation/` reproducing PyLCM-vs-Fortran/SAM plots + the `E_S09` explanation | `tests/`, `validation/` |
| **D. CI** | `.github/workflows/ci.yml` — install env, run pytest on push/PR (py3.11 + 3.12 matrix) | `.github/` |
| **E. Ensemble + speed** | Promote `ensemble_comparison.py` → `ensemble.py` module/CLI: N members, varying RNG seed, parallel across cores (joblib), mean+quantile envelope plots. Profile single run, expand `@njit` on hot spots, document measured before/after | `ensemble.py`, hot modules |
| **F. Placeholder gate** | Hard-gate `entrainment.py` behind explicit `experimental=True` + loud warning + docstring "not physically validated" | `entrainment.py`, README |
| **G. README** (last) | Rewrite: badges, one-command install, quickstart, the two notebooks, scientific basis + refs, validation summary, known limitations, how-to-cite, contact | `README.md` |

## 6. Testing strategy

Test **physics invariants**, not golden numbers (invariants survive Monte Carlo variance):

- **Collision conserves** total water mass and total number·mass before/after a collision step.
- **Köhler** equilibrium radius is positive and behaves monotonically with the expected inputs.
- **Condensation** produces no NaN/Inf over a short integration from a sane initial state.
- **`aero_init` unit conversion** sanity: `mu = log(mu_um·1e-6)`, `sigma = log(sigma_geo)`
  (guards the documented overflow-on-raw-values crash).
- **Ensemble reproducibility**: same seed → identical output; different seeds → bounded spread.

CI runs the suite on every push/PR so stability is provable going forward.

## 7. Dependency order

`Setup (commit + worktree)`
→ `A Hygiene` + `B Packaging` (early, so junk does not return)
→ parallel `C Tests` + `E Ensemble/speed` + `F Gate`
→ `D CI` (needs tests + env)
→ `G README` (needs CI badge + final commands)
→ merge `release/v1.0` → `main` + tag `v1.0.0`.

## 8. Decisions (resolved)

- **License:** MIT.
- **Worktree path:** `../PyLCM_v1.0`.
- **Final merge:** merge `release/v1.0` → `main` and tag `v1.0.0`.
- **Entrainment:** deferred to v1.1; hard-gated as EXPERIMENTAL in v1.0.
- **Performance scope:** profile + expand JIT + parallel ensemble; no risky rewrites.

## 9. Acceptance criteria

- Fresh clone + one documented command yields a working environment.
- `pytest` passes locally and in CI; CI badge green in README.
- Collision-validation notebook runs end-to-end and renders the comparison + conclusion.
- `ensemble.py` runs an N-member ensemble in parallel with a measured wall-clock win;
  before/after numbers recorded.
- `entrainment` cannot be invoked without an explicit experimental opt-in.
- No build artifacts tracked; `LICENSE` present; `main` tagged `v1.0.0`.
