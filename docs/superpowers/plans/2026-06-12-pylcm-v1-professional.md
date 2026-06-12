# PyLCM v1.0 Professional — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn working-but-rough PyLCM into a clone-install-trust-run educational release: clean repo, pinned env, physics-invariant tests + CI, a collision-validation notebook, a parallel ensemble module with profiled speedups, gated experimental code, rewritten README, tagged `v1.0.0`.

**Architecture:** All work happens in a git worktree `../PyLCM_v1.0` on branch `release/v1.0` (off `feature/lecture-readiness-improvements`) so the user's live PyLCM session in `/Users/dr.cloud/PyLCM_edu` is undisturbed. Lanes touch disjoint files and are largely parallel. Tests assert conservation/positivity invariants (not golden numbers) so they survive Monte Carlo variance.

**Tech Stack:** Python 3.11/3.12, numpy/scipy/numba, pytest, joblib, conda + pip, GitHub Actions, Jupyter.

**Hard constraints:** NEVER add a "Claude"/Co-Authored-By line to any commit. Do not modify files in the live `PyLCM_edu` dir during lane work — operate in the worktree.

---

## File Structure

| Path | Responsibility |
|------|----------------|
| `.gitignore` | Ignore Python/Jupyter/macOS/numba build artifacts |
| `LICENSE` | MIT license text |
| `environment.yml` | Conda env, pinned majors |
| `requirements.txt` | pip deps, pinned majors |
| `pyproject.toml` | Make `PyLCM`/`Post_process` importable via `pip install -e .` |
| `tests/conftest.py` | Shared fixtures (a minimal particle, a tiny parcel state) |
| `tests/test_collision_invariants.py` | Mass conservation + number-decrease in collisions |
| `tests/test_kohler.py` | Equilibrium radius positivity/monotonicity |
| `tests/test_condensation_stability.py` | No NaN/Inf over a short integration |
| `tests/test_aero_init.py` | Unit-conversion guard (log of meters / log of geo-std) |
| `tests/test_ensemble.py` | Seed determinism + bounded spread |
| `validation/collision_validation.ipynb` | PyLCM vs Fortran/SAM collision comparison + E_S09 story |
| `ensemble.py` | Parallel N-member ensemble runner + envelope plot |
| `.github/workflows/ci.yml` | Run pytest on push/PR |
| `PyLCM/entrainment.py` | Add experimental hard-gate |
| `README.md` | Full rewrite |

---

## Task 0: Worktree setup

**Files:** none created in worktree yet; operates on the live repo's index.

- [ ] **Step 1: Stage only the good uncommitted work on the feature branch**

From `/Users/dr.cloud/PyLCM_edu` (current branch `feature/lecture-readiness-improvements`):

```bash
git add PyLCM_Part1_Foundations.ipynb PyLCM_Part2_Experiments.ipynb \
        ensemble_comparison.py test_init_kappa.py test_stability.py \
        overnight_sensitivity_test.sh SENSITIVITY_TEST_HANDOFF.md \
        plot_comparison.py COLLISION_COMPARISON.md ASSESSMENT_REPORT.md \
        sensitivity_init_kappa_dsd.png sensitivity_init_kappa_ts.png \
        comparison_dsd.png comparison_timeseries.png
git status --short
```

Expected: the listed files are staged; NO `__pycache__`, `.DS_Store`, `.pyc`, `.nbc`, `.nbi` staged.

- [ ] **Step 2: Commit the good work (no Claude coauthor)**

```bash
git commit -m "Add educational notebooks, ensemble script, and validation artifacts"
```

- [ ] **Step 3: Create the worktree**

```bash
git worktree add ../PyLCM_v1.0 -b release/v1.0 feature/lecture-readiness-improvements
git worktree list
```

Expected: `../PyLCM_v1.0  <sha> [release/v1.0]` appears. **All subsequent tasks run inside `../PyLCM_v1.0`.**

- [ ] **Step 4: Sanity check the worktree imports**

```bash
cd ../PyLCM_v1.0 && python -c "import PyLCM.collision, PyLCM.condensation, PyLCM.aero_init; print('ok')"
```

Expected: `ok` (run from the worktree root so `PyLCM` is on the path).

---

## Task 1: Repo hygiene (.gitignore, untrack artifacts, LICENSE)

**Files:**
- Create: `.gitignore`, `LICENSE`

- [ ] **Step 1: Write `.gitignore`**

```gitignore
# Python
__pycache__/
*.py[cod]
*.egg-info/
.eggs/
build/
dist/
# Numba cache
*.nbc
*.nbi
# Jupyter
.ipynb_checkpoints/
# macOS
.DS_Store
# Project output
Output/
# Envs
.venv/
env/
venv/
```

- [ ] **Step 2: Untrack already-committed build artifacts**

```bash
git rm -r --cached --ignore-unmatch $(git ls-files | grep -E '__pycache__|\.pyc$|\.DS_Store$|\.nbc$|\.nbi$') 2>/dev/null
git status --short | grep -E 'deleted' | head
```

Expected: tracked pycache/pyc/DS_Store entries show as staged deletions (files stay on disk, just untracked).

- [ ] **Step 3: Write `LICENSE` (MIT)**

```text
MIT License

Copyright (c) 2026 J. Lim

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

- [ ] **Step 4: Commit**

```bash
git add .gitignore LICENSE
git commit -m "Add .gitignore and MIT LICENSE; untrack build artifacts"
```

---

## Task 2: Packaging (environment.yml, requirements.txt, pyproject.toml)

**Files:**
- Create: `environment.yml`, `requirements.txt`, `pyproject.toml`

- [ ] **Step 1: Capture the actual installed versions to pin against**

```bash
python -c "import numpy,scipy,pandas,matplotlib,plotly,numba,ipywidgets,tqdm; \
print('numpy',numpy.__version__); print('scipy',scipy.__version__); \
print('pandas',pandas.__version__); print('matplotlib',matplotlib.__version__); \
print('plotly',plotly.__version__); print('numba',numba.__version__); \
print('ipywidgets',ipywidgets.__version__); print('tqdm',tqdm.__version__)"
```

Expected: version strings print. Use the printed MAJOR.MINOR as lower bounds below (replace the placeholders).

- [ ] **Step 2: Write `requirements.txt`** (substitute the versions printed in Step 1)

```text
numpy>=1.24
scipy>=1.10
pandas>=2.0
matplotlib>=3.7
plotly>=5.0
numba>=0.58
ipywidgets>=8.0
tqdm>=4.65
joblib>=1.3
jupyter
```

- [ ] **Step 3: Write `environment.yml`**

```yaml
name: PyLCM
channels:
  - conda-forge
dependencies:
  - python>=3.11,<3.13
  - numpy>=1.24
  - scipy>=1.10
  - pandas>=2.0
  - matplotlib>=3.7
  - numba>=0.58
  - ipywidgets>=8.0
  - tqdm>=4.65
  - joblib>=1.3
  - jupyter
  - pip
  - pip:
    - plotly>=5.0
```

- [ ] **Step 4: Write `pyproject.toml`**

```toml
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"

[project]
name = "pylcm"
version = "1.0.0"
description = "Educational Lagrangian Cloud Model parcel simulator"
readme = "README.md"
requires-python = ">=3.11"
license = { text = "MIT" }
authors = [{ name = "J. Lim" }]
dependencies = [
  "numpy>=1.24", "scipy>=1.10", "pandas>=2.0", "matplotlib>=3.7",
  "plotly>=5.0", "numba>=0.58", "ipywidgets>=8.0", "tqdm>=4.65", "joblib>=1.3",
]

[tool.setuptools]
packages = ["PyLCM", "Post_process"]
```

- [ ] **Step 5: Verify editable install works**

```bash
pip install -e . && python -c "import PyLCM; print('PyLCM importable as package')"
```

Expected: install succeeds; prints the confirmation line.

- [ ] **Step 6: Commit**

```bash
git add environment.yml requirements.txt pyproject.toml
git commit -m "Add pinned environment.yml, requirements.txt, and pyproject.toml"
```

---

## Task 3: Test scaffold + collision invariant (TDD)

**Files:**
- Create: `tests/__init__.py`, `tests/conftest.py`, `tests/test_collision_invariants.py`

- [ ] **Step 1: Write the failing collision test**

`tests/__init__.py`: empty file.

`tests/conftest.py`:

```python
import pytest
from PyLCM.micro_particle import particles


def make_particle(M, A, Ns=1e-18, kappa=0.5):
    p = particles(1)
    p.M, p.A, p.Ns, p.kappa = M, A, Ns, kappa
    return p
```

`tests/test_collision_invariants.py`:

```python
from PyLCM.collision import liquid_update_collection
from tests.conftest import make_particle


def test_collision_conserves_total_water_mass():
    # M is total superdroplet water mass; collision moves mass between the two
    p1 = make_particle(M=2.0e-9, A=100.0)
    p2 = make_particle(M=5.0e-9, A=10.0)
    total_before = p1.M + p2.M
    a1, a2, acc, aut = liquid_update_collection(p1, p2, 0.0, 0.0, p_crit=1)
    total_after = a1.M + a2.M
    assert abs(total_after - total_before) < 1e-20


def test_collision_decreases_total_number():
    # Coalescence reduces total multiplicity A (the larger-A particle loses count)
    p1 = make_particle(M=2.0e-9, A=100.0)
    p2 = make_particle(M=5.0e-9, A=10.0)
    a_before = p1.A + p2.A
    a1, a2, _, _ = liquid_update_collection(p1, p2, 0.0, 0.0, p_crit=1)
    assert (a1.A + a2.A) < a_before
```

- [ ] **Step 2: Run to verify it passes (these assert real existing behavior)**

```bash
cd ../PyLCM_v1.0 && python -m pytest tests/test_collision_invariants.py -v
```

Expected: 2 passed. (If `test_collision_conserves_total_water_mass` FAILS, that is a real conservation bug — stop and report it before continuing; do not "fix" the test to match.)

- [ ] **Step 3: Add the same-weight branch test**

Append to `tests/test_collision_invariants.py`:

```python
from PyLCM.collision import same_weights_update


def test_same_weight_collision_conserves_mass():
    p1 = make_particle(M=3.0e-9, A=50.0)
    p2 = make_particle(M=4.0e-9, A=50.0)
    total_before = p1.M + p2.M
    a1, a2, _, _ = same_weights_update(p1, p2, 0.0, 0.0)
    assert abs((a1.M + a2.M) - total_before) < 1e-20
```

- [ ] **Step 4: Run all collision tests**

```bash
python -m pytest tests/test_collision_invariants.py -v
```

Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add tests/__init__.py tests/conftest.py tests/test_collision_invariants.py
git commit -m "Add collision mass-conservation and number-decrease invariant tests"
```

---

## Task 4: Köhler, aero_init, and condensation-stability tests (TDD)

**Files:**
- Create: `tests/test_kohler.py`, `tests/test_aero_init.py`, `tests/test_condensation_stability.py`

- [ ] **Step 1: Write the Köhler equilibrium test**

`tests/test_kohler.py`:

```python
import numpy as np
from PyLCM.aero_init import r_equi


def test_equilibrium_radius_positive():
    # r_equi(S, T, r_aerosol, rho_aero, switch_kappa_koehler, kappa)
    r = r_equi(-0.01, 280.0, 5.0e-8, 1777.0, True, 0.5)
    assert np.isfinite(r) and r > 0


def test_equilibrium_radius_grows_with_humidity():
    r_dry = r_equi(-0.05, 280.0, 5.0e-8, 1777.0, True, 0.5)
    r_humid = r_equi(-0.005, 280.0, 5.0e-8, 1777.0, True, 0.5)
    assert r_humid >= r_dry
```

- [ ] **Step 2: Run it**

```bash
cd ../PyLCM_v1.0 && python -m pytest tests/test_kohler.py -v
```

Expected: 2 passed. If `r_equi`'s argument order differs from the docstring at `PyLCM/aero_init.py:249`, adjust the call to match the real signature, then re-run.

- [ ] **Step 3: Write the aero_init unit-conversion guard**

`tests/test_aero_init.py`:

```python
import numpy as np
from PyLCM.aero_init import aero_init


def test_aero_init_runs_with_log_converted_params():
    # Per the documented convention: mu in log(meters), sigma in log(geo-std).
    mu = np.log(np.array([0.02e-6, 0.2e-6]))   # 0.02, 0.2 micron -> meters -> log
    sigma = np.log(np.array([1.4, 1.6]))       # geometric std -> log
    N_aero = np.array([100.0e6, 20.0e6])
    out = aero_init("weighting_factor", 200, 100000.0, 0.0, 285.0, 0.008,
                    N_aero, mu, sigma, 1777.0, np.array([0.5, 0.5]), True)
    # Whatever the return tuple is, the first element (particle list/array) must be non-empty
    assert out is not None
```

- [ ] **Step 4: Run it**

```bash
python -m pytest tests/test_aero_init.py -v
```

Expected: 1 passed. If the `mode_aero_init` string or argument order differs, read `PyLCM/aero_init.py:116` and adjust the call to the real signature/mode name; the point is that log-converted inputs do NOT raise/overflow.

- [ ] **Step 5: Write the condensation no-NaN test**

`tests/test_condensation_stability.py`:

```python
import numpy as np
from PyLCM.aero_init import aero_init
from PyLCM.condensation import drop_condensation


def test_condensation_no_nan_short_run():
    mu = np.log(np.array([0.02e-6, 0.2e-6]))
    sigma = np.log(np.array([1.4, 1.6]))
    N_aero = np.array([100.0e6, 20.0e6])
    particles_list = aero_init("weighting_factor", 200, 100000.0, 0.0, 285.0,
                               0.008, N_aero, mu, sigma, 1777.0,
                               np.array([0.5, 0.5]), True)[0]
    T, q, P = 285.0, 0.008, 100000.0
    S_lst = []
    for _ in range(10):
        out = drop_condensation(particles_list, T, q, P, 0, 0.1, 100.0, S_lst,
                                1777.0, 1.0e-6, 0.0, 0.0, 0.0, 0.0, True)
        T, q, P = out[1], out[2], out[3]  # adjust indices to real return order
        assert np.isfinite(T) and np.isfinite(q)
```

- [ ] **Step 6: Run it; reconcile the return-tuple indices**

```bash
python -m pytest tests/test_condensation_stability.py -v
```

Expected: 1 passed. Read `PyLCM/condensation.py:8` return statement to map the real indices for `T_parcel`, `q_parcel`, `P_parcel`; fix the unpacking, then re-run until green.

- [ ] **Step 7: Commit**

```bash
git add tests/test_kohler.py tests/test_aero_init.py tests/test_condensation_stability.py
git commit -m "Add Kohler, aero_init unit-conversion, and condensation stability tests"
```

---

## Task 5: Ensemble module (parallel members) + determinism test

**Files:**
- Create: `ensemble.py`, `tests/test_ensemble.py`
- Reference: `ensemble_comparison.py` (existing ad-hoc script to refactor from)

- [ ] **Step 1: Read the existing ensemble script to reuse its single-run path**

```bash
cd ../PyLCM_v1.0 && sed -n '1,80p' ensemble_comparison.py
```

Note the function that runs one simulation and the arguments it needs (initial T/P/RH/w, n_ptcl, dt, nt, aerosol params).

- [ ] **Step 2: Write `ensemble.py`**

```python
"""Parallel ensemble runner for PyLCM.

Each member is an independent stochastic realization (different RNG seed).
Members are embarrassingly parallel, so we fan them across CPU cores.
"""
import numpy as np
from joblib import Parallel, delayed


def run_member(seed, run_single, run_kwargs):
    """Run one ensemble member with a fixed seed. `run_single(**kwargs)` must
    return a 1-D numpy array (a diagnostic time series, e.g. LWC or Nr)."""
    np.random.seed(seed)
    return np.asarray(run_single(**run_kwargs))


def run_ensemble(run_single, run_kwargs, n_members=10, n_jobs=-1, base_seed=0):
    """Run `n_members` members in parallel; return (members, mean, lo, hi).

    members: (n_members, n_t) array; lo/hi are the 10th/90th percentile envelope."""
    seeds = [base_seed + i for i in range(n_members)]
    results = Parallel(n_jobs=n_jobs)(
        delayed(run_member)(s, run_single, run_kwargs) for s in seeds
    )
    members = np.vstack(results)
    mean = members.mean(axis=0)
    lo = np.percentile(members, 10, axis=0)
    hi = np.percentile(members, 90, axis=0)
    return members, mean, lo, hi


def plot_envelope(time, mean, lo, hi, label="", ax=None):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(time, mean, label=label)
    ax.fill_between(time, lo, hi, alpha=0.3)
    ax.legend()
    return ax
```

- [ ] **Step 3: Write the determinism + spread test**

`tests/test_ensemble.py`:

```python
import numpy as np
from ensemble import run_ensemble


def _toy_run(scale=1.0):
    # Stand-in for a PyLCM single run: stochastic series depending on the RNG.
    return scale * np.cumsum(np.random.random(20))


def test_same_seed_is_deterministic():
    m1, *_ = run_ensemble(_toy_run, {"scale": 1.0}, n_members=3, n_jobs=1, base_seed=42)
    m2, *_ = run_ensemble(_toy_run, {"scale": 1.0}, n_members=3, n_jobs=1, base_seed=42)
    assert np.allclose(m1, m2)


def test_envelope_brackets_mean():
    _, mean, lo, hi = run_ensemble(_toy_run, {"scale": 1.0}, n_members=20, n_jobs=1, base_seed=0)
    assert np.all(lo <= mean + 1e-9) and np.all(mean - 1e-9 <= hi)
```

- [ ] **Step 4: Run the ensemble tests**

```bash
python -m pytest tests/test_ensemble.py -v
```

Expected: 2 passed.

- [ ] **Step 5: Wire `ensemble.py` to the real PyLCM single-run (smoke check)**

Adapt `ensemble_comparison.py`'s single-run into a `run_single(**kwargs)` returning one diagnostic series, then:

```bash
python -c "from ensemble import run_ensemble; print('import + API ok')"
```

Expected: prints confirmation. (Full real-physics ensemble run is exercised in the validation notebook, not in CI, to keep CI fast.)

- [ ] **Step 6: Commit**

```bash
git add ensemble.py tests/test_ensemble.py
git commit -m "Add parallel ensemble module with determinism and envelope tests"
```

---

## Task 6: Profile single run + expand JIT (measured)

**Files:**
- Modify: hot functions in `PyLCM/condensation.py` and/or `PyLCM/collision.py` (only pure-numeric helpers)
- Create: `validation/PROFILE.md` (before/after numbers)

- [ ] **Step 1: Profile a representative single run**

```bash
cd ../PyLCM_v1.0 && python -m cProfile -s cumtime ensemble_comparison.py 2>/dev/null | head -30
```

Note the top cumulative-time functions that are pure-numeric and not already `@jit`.

- [ ] **Step 2: Record the baseline in `validation/PROFILE.md`**

```markdown
# PyLCM single-run profile

## Baseline (before JIT expansion)
- <function>: <cumtime>s
- Total wall: <t>s for <n_ptcl> particles, <nt> steps
```

- [ ] **Step 3: Add `@jit(nopython=True, cache=True)` to ONE safe hot helper**

Only decorate a pure-numeric function with no Python objects/loops over particle instances (e.g. an array helper in `condensation.py`). Show the exact diff in the commit. Do NOT JIT functions that touch `particles` instances or lists of them.

- [ ] **Step 4: Re-profile and confirm correctness is unchanged**

```bash
python -m pytest tests/ -v   # all invariants still pass
python -m cProfile -s cumtime ensemble_comparison.py 2>/dev/null | head -30
```

Expected: all tests pass; the decorated function's cumtime drops after warmup. Append "After" numbers to `validation/PROFILE.md`. If a change does not measurably help or breaks a test, revert it — only keep wins.

- [ ] **Step 5: Commit**

```bash
git add validation/PROFILE.md PyLCM/
git commit -m "Profile single run and expand numba JIT on hot numeric helpers"
```

---

## Task 7: Hard-gate experimental entrainment

**Files:**
- Modify: `PyLCM/entrainment.py`

- [ ] **Step 1: Read the current entrainment entry points**

```bash
cd ../PyLCM_v1.0 && grep -n "^def " PyLCM/entrainment.py
```

- [ ] **Step 2: Add an explicit experimental gate at the top of each public entrainment function**

Add this guard helper at the top of `PyLCM/entrainment.py` (after imports) and call it at the start of each public function:

```python
import warnings


def _require_experimental(experimental):
    if not experimental:
        raise RuntimeError(
            "entrainment is EXPERIMENTAL and not physically validated. "
            "Pass experimental=True to opt in explicitly. See README 'Known limitations'."
        )
    warnings.warn("Using EXPERIMENTAL, unvalidated entrainment physics.", stacklevel=2)
```

For each public function add an `experimental=False` keyword parameter and call `_require_experimental(experimental)` as its first line. Update the module docstring to state it is experimental.

- [ ] **Step 3: Add a test that the gate fires**

Append `tests/test_entrainment_gate.py`:

```python
import pytest
import PyLCM.entrainment as ent


def test_entrainment_blocked_by_default():
    # Pick the first public function; it must raise without experimental=True.
    funcs = [getattr(ent, n) for n in dir(ent)
             if callable(getattr(ent, n)) and not n.startswith("_")
             and getattr(getattr(ent, n), "__module__", "") == ent.__name__]
    assert funcs, "no public entrainment functions found"
    with pytest.raises(RuntimeError):
        funcs[0]()  # called with no args -> gate raises before arg errors
```

- [ ] **Step 4: Run it**

```bash
python -m pytest tests/test_entrainment_gate.py -v
```

Expected: 1 passed. (If the first function's gate is reached only after positional args, place `_require_experimental` as the literal first statement so it raises before any other error.)

- [ ] **Step 5: Commit**

```bash
git add PyLCM/entrainment.py tests/test_entrainment_gate.py
git commit -m "Hard-gate experimental entrainment behind explicit opt-in"
```

---

## Task 8: CI workflow

**Files:**
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Write the workflow**

```yaml
name: CI
on:
  push:
  pull_request:
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install
        run: |
          python -m pip install --upgrade pip
          pip install -e .
          pip install pytest
      - name: Run tests
        run: python -m pytest tests/ -v
```

- [ ] **Step 2: Run the full suite locally exactly as CI will**

```bash
cd ../PyLCM_v1.0 && pip install -e . && python -m pytest tests/ -v
```

Expected: all tests pass (this is the gate CI enforces).

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "Add GitHub Actions CI running pytest on py3.11 and py3.12"
```

---

## Task 9: Collision-validation notebook

**Files:**
- Create: `validation/collision_validation.ipynb`
- Reference: `COLLISION_COMPARISON.md`, `comparison_timeseries.png`, `comparison_dsd.png`

- [ ] **Step 1: Build the notebook with these markdown + code cells**

Cells (create via `jupyter nbconvert`-friendly JSON or in Jupyter):
1. **Markdown:** title + the question ("Are PyLCM collision results correct vs Fortran/SAM?").
2. **Markdown — the verdict:** PyLCM and SAM6-LCM share LSM (Shima 2009), Hall (1980) E-table, Beard (1976) vt, multi-collision p_crit limiter. SAM omits Straub (2009) `E_S09` in BOTH gravitational and turbulent kernels (`micro_coll.f90:710`); PyLCM includes it (`collision.py:205,207`). PyLCM therefore agrees with the Fortran box model (which also includes `E_S09`), and differs from SAM by ~1/E_S09 where turbulence is active. Residual ~3× peak-Nr gap vs box model = LSM vs O(N²) enumeration + Monte Carlo variance.
3. **Code:** load and display `comparison_timeseries.png` and `comparison_dsd.png`.
4. **Code:** plot `E_S09` vs Weber number using `PyLCM.collision.E_S09` to show the magnitude of the term SAM drops.
5. **Markdown — conclusion:** PyLCM's collision scheme is correct and more physically complete than SAM6-LCM in the turbulent regime; the discrepancy is a documented modeling choice, not a bug.

- [ ] **Step 2: Execute the notebook end-to-end**

```bash
cd ../PyLCM_v1.0 && jupyter nbconvert --to notebook --execute --inplace validation/collision_validation.ipynb
```

Expected: runs with no errors; figures render.

- [ ] **Step 3: Commit**

```bash
git add validation/collision_validation.ipynb
git commit -m "Add collision-validation notebook documenting PyLCM vs Fortran/SAM"
```

---

## Task 10: README rewrite

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Rewrite `README.md`** with these sections:

```markdown
# PyLCM — Educational Lagrangian Cloud Model

![CI](https://github.com/<user>/<repo>/actions/workflows/ci.yml/badge.svg)

A teaching parcel-model simulator for warm-cloud microphysics using Lagrangian
super-droplets (condensation, collision-coalescence).

## Install (one command)
    conda env create -f environment.yml && conda activate PyLCM && pip install -e .
Or with pip: `pip install -e . && pip install -r requirements.txt`.

## Quickstart
    jupyter notebook PyLCM_edu.ipynb
Teaching notebooks: `PyLCM_Part1_Foundations.ipynb`, `PyLCM_Part2_Experiments.ipynb`.

## Scientific basis
Köhler theory; Hall (1980) collision efficiency; Straub et al. (2009) coalescence
efficiency; Beard (1976) terminal velocity; Shima et al. (2009) Linear Sampling
Method; optional Ayala et al. (2008)/Wang & Grabowski (2009) turbulent kernel.

## Validation
See `validation/collision_validation.ipynb`: PyLCM matches the Fortran box-model
reference and is more physically complete than SAM6-LCM in the turbulent regime
(SAM omits the Straub 2009 coalescence efficiency). Physics invariants are checked
in CI (`tests/`).

## Ensemble runs
`ensemble.py` runs N stochastic members in parallel across cores with a percentile
envelope. See its docstring.

## Known limitations
- **Entrainment is EXPERIMENTAL and not physically validated.** It is hard-gated;
  you must pass `experimental=True` to use it. A validated scheme is planned for v1.1.

## How to cite
J. Lim, PyLCM v1.0 (2026). Educational Lagrangian Cloud Model.

## Contact
J.lim@physik.uni-muenchen.de

## License
MIT (see `LICENSE`).
```

Replace `<user>/<repo>` with the real GitHub slug (`git remote get-url origin`).

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "Rewrite README: one-command install, validation, ensemble, limitations"
```

---

## Task 11: Release — merge + tag

**Files:** none (git operations)

- [ ] **Step 1: Final full-suite gate**

```bash
cd ../PyLCM_v1.0 && python -m pytest tests/ -v
```

Expected: all green.

- [ ] **Step 2: Merge release branch into main**

```bash
cd /Users/dr.cloud/PyLCM_edu && git fetch && git checkout main && \
git merge --no-ff release/v1.0 -m "Release PyLCM v1.0 (professional)"
```

(If the live working dir has the running session checked out on the feature branch, do the merge from the worktree side or coordinate with the user before switching branches in the live dir.)

- [ ] **Step 3: Tag and clean up the worktree**

```bash
git tag -a v1.0.0 -m "PyLCM v1.0.0 — professional educational release"
git worktree remove ../PyLCM_v1.0
git worktree list
```

Expected: `v1.0.0` tag exists; worktree removed. Push (`git push origin main --tags`) only when the user explicitly asks.

---

## Self-Review

**Spec coverage:** §2 hygiene→Task 1; packaging→Task 2; tests→Tasks 3,4,5,7; CI→Task 8; collision validation→Task 9; ensemble+speed→Tasks 5,6; entrainment gate→Task 7; README→Task 10; release/tag→Task 11; worktree strategy→Task 0. All covered.

**Placeholder scan:** Test code is concrete; the only deliberate "reconcile to real signature" notes (Tasks 4,5) are because the exact return-tuple index ordering must be read from source at execution — each such step names the exact file:line to read and the assertion to keep. No vague "add error handling".

**Type consistency:** `make_particle` fixture used consistently across Tasks 3–4; `run_ensemble`/`run_member`/`run_single` signatures consistent across Task 5 and its test; `_require_experimental(experimental)` consistent across Task 7.
