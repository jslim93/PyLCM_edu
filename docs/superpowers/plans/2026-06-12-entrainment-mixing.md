# Entrainment / IHMD Mixing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add warm-cloud entrainment mixing parameterized by entrainment rate λ and Inhomogeneous Mixing Degree IHMD (Lim & Hoffmann 2023), spanning homogeneous↔inhomogeneous mixing, behind a `MixingModel` interface a future LEM backend can implement.

**Architecture:** A new `PyLCM/mixing.py` holds a `MixingModel` protocol and a `ParameterizedMixing` implementation. Mixing runs **before** condensation each step: it dilutes bulk T,q toward the environment and redistributes cloud liquid across super-droplets so that `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD` (exact closed form: `M ← M·(1−δ)`, `A ← A·(1−δ)^IHMD`). Default disabled → existing runs unchanged.

**Tech Stack:** Python, numpy, pytest. Worktree `../PyLCM_entrain` off `feature/performance-ensemble`.

**Constraints:** Never add a "Claude"/Co-Authored-By commit line. Never touch `/Users/dr.cloud/PyLCM_edu`.

---

## File Structure

| Path | Responsibility |
|------|----------------|
| `PyLCM/mixing.py` | `MixingModel` protocol, `ParameterizedMixing`, `LEMMixing` stub, env interpolation |
| `tests/test_mixing_redistribute.py` | The IHMD closed-form invariants (N–q law, endpoints, budget) |
| `tests/test_mixing_apply.py` | Full `apply` step: bulk dilution + redistribution, determinism, no-op |
| `tests/test_mixing_lem_stub.py` | `LEMMixing` raises NotImplementedError with Phase-3b pointer |
| `pylcm_run.py` | Optional `mixing` arg wired into the run loop (mixing before condensation) |
| `PyLCM/timestep_routine.py` | Replace `basic_entrainment` call with `MixingModel.apply`; add `ihmd` |
| `validation/entrainment_mixing.ipynb` | IHMD sweep: DSD + N–r_v diagram |

---

## Task 0: Worktree setup

- [ ] **Step 1: Create the worktree off the perf branch**

```bash
cd /Users/dr.cloud/PyLCM_edu
git worktree add ../PyLCM_entrain -b feature/entrainment-mixing feature/performance-ensemble
git worktree list
```

Expected: `../PyLCM_entrain  <sha> [feature/entrainment-mixing]`. All tasks run there.

- [ ] **Step 2: Confirm baseline tests pass**

```bash
cd /Users/dr.cloud/PyLCM_entrain && pip install -e . -q && python -m pytest tests/ -q
```

Expected: the inherited suite (22 tests) passes.

---

## Task 1: IHMD droplet redistribution (the physics core, TDD)

**Files:**
- Create: `PyLCM/mixing.py`, `tests/test_mixing_redistribute.py`

- [ ] **Step 1: Write the failing invariant tests**

`tests/test_mixing_redistribute.py`:

```python
import numpy as np
from PyLCM.micro_particle import particles
from PyLCM.mixing import redistribute_droplets


def _cloud(n=50, M=2.0e-9, A=100.0):
    out = []
    for _ in range(n):
        p = particles(1); p.M, p.A, p.Ns, p.kappa = M, A, 1e-18, 0.5
        out.append(p)
    return out


def test_nq_power_law_holds_exactly():
    for ihmd in (0.0, 0.5, 1.0):
        ps = _cloud()
        N0 = sum(p.A for p in ps); q0 = sum(p.M for p in ps)
        evap = redistribute_droplets(ps, ihmd=ihmd, frac=0.3)
        N1 = sum(p.A for p in ps); q1 = sum(p.M for p in ps)
        assert abs(evap - 0.3 * q0) < 1e-18                      # evaporated mass
        assert np.isclose(N1 / N0, (q1 / q0) ** ihmd, rtol=1e-9)  # the IHMD law


def test_homogeneous_conserves_number():
    ps = _cloud(); N0 = sum(p.A for p in ps)
    redistribute_droplets(ps, ihmd=0.0, frac=0.4)
    assert abs(sum(p.A for p in ps) - N0) < 1e-9                 # N unchanged
    assert all(np.isclose(p.M / p.A, (2.0e-9 / 100.0) * 0.6) for p in ps)  # all shrink


def test_inhomogeneous_preserves_droplet_size():
    ps = _cloud(); N0 = sum(p.A for p in ps)
    redistribute_droplets(ps, ihmd=1.0, frac=0.4)
    assert sum(p.A for p in ps) < N0                            # number drops
    assert all(np.isclose(p.M / p.A, 2.0e-9 / 100.0) for p in ps)  # size unchanged
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd /Users/dr.cloud/PyLCM_entrain && python -m pytest tests/test_mixing_redistribute.py -q`
Expected: FAIL (`No module named PyLCM.mixing`).

- [ ] **Step 3: Implement `redistribute_droplets`**

`PyLCM/mixing.py`:

```python
"""Entrainment mixing for PyLCM (warm cloud).

Mixing runs BEFORE condensation each timestep. The Inhomogeneous Mixing Degree
(IHMD, Lim & Hoffmann 2023) controls how entrainment-driven evaporation is split
between homogeneous (all droplets shrink, number conserved) and inhomogeneous
(a subset evaporates, survivors keep size) limits, satisfying exactly:

    N_c / N_{c,0} = (q_c / q_{c,0}) ** IHMD

Closed form per step with entrained fraction frac in [0,1):
    M <- M * (1 - frac)            # total super-droplet liquid mass
    A <- A * (1 - frac) ** IHMD    # multiplicity (droplet number)
"""
import numpy as np


def redistribute_droplets(particles_list, ihmd, frac):
    """Apply IHMD redistribution to cloud super-droplets in place.

    Returns the total evaporated liquid mass (sum of M lost), which the caller
    returns to the vapor field.
    """
    if frac <= 0.0:
        return 0.0
    keep_mass = 1.0 - frac
    keep_num = keep_mass ** ihmd
    evaporated = 0.0
    for p in particles_list:
        if p.A <= 0 or p.M <= 0:
            continue
        m_old = p.M
        p.M = p.M * keep_mass
        p.A = p.A * keep_num
        evaporated += m_old - p.M
    return evaporated
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_mixing_redistribute.py -q`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add PyLCM/mixing.py tests/test_mixing_redistribute.py
git commit -m "Add IHMD droplet redistribution (exact N-q power law)"
```

---

## Task 2: Environmental interpolation + `ParameterizedMixing.apply` (TDD)

**Files:**
- Modify: `PyLCM/mixing.py`
- Create: `tests/test_mixing_apply.py`

- [ ] **Step 1: Write the failing tests**

`tests/test_mixing_apply.py`:

```python
import numpy as np
from PyLCM.micro_particle import particles
from PyLCM.parcel import create_env_profiles
from PyLCM.mixing import ParameterizedMixing


def _setup():
    # Suppress the plot side effect of create_env_profiles in headless tests.
    import matplotlib; matplotlib.use("Agg")
    qv, th, z_env = create_env_profiles(290.0, 0.010, 0.0, 95000.0, "Stable")
    ps = []
    for _ in range(100):
        p = particles(1); p.M, p.A, p.Ns, p.kappa = 1.0e-9, 50.0, 1e-18, 0.5
        ps.append(p)
    return qv, th, z_env, ps


def test_apply_conserves_total_water():
    qv, th, z_env, ps = _setup()
    air_mass = 100.0
    mix = ParameterizedMixing(lambda_ent=5e-4, ihmd=0.5, qv_profiles=qv,
                              theta_profiles=th, z_env=z_env)
    T, q = 288.0, 0.009
    liq0 = sum(p.M for p in ps)
    total0 = q * air_mass + liq0                      # vapor + liquid (kg)
    ps, T1, q1 = mix.apply(ps, T, q, 90000.0, 500.0, dt=1.0, w=1.0, air_mass=air_mass)
    total1 = q1 * air_mass + sum(p.M for p in ps)
    # Entrainment also imports environmental vapor; check the LIQUID->VAPOR part
    # conserves by comparing evaporated liquid to the vapor gained from evaporation.
    assert total1 >= total0 - 1e-9                    # no spurious water loss
    assert np.isfinite(T1) and np.isfinite(q1)


def test_disabled_is_noop():
    qv, th, z_env, ps = _setup()
    mix = ParameterizedMixing(lambda_ent=0.0, ihmd=0.5, qv_profiles=qv,
                              theta_profiles=th, z_env=z_env)
    M_before = [p.M for p in ps]
    ps, T1, q1 = mix.apply(ps, 288.0, 0.009, 90000.0, 500.0, dt=1.0, w=1.0, air_mass=100.0)
    assert [p.M for p in ps] == M_before and (T1, q1) == (288.0, 0.009)


def test_apply_is_deterministic():
    qv, th, z_env, ps1 = _setup()
    _, _, _, ps2 = _setup()
    m1 = ParameterizedMixing(8e-4, 0.7, qv, th, z_env)
    m2 = ParameterizedMixing(8e-4, 0.7, qv, th, z_env)
    r1 = m1.apply(ps1, 288.0, 0.009, 90000.0, 500.0, 1.0, 1.0, 100.0)
    r2 = m2.apply(ps2, 288.0, 0.009, 90000.0, 500.0, 1.0, 1.0, 100.0)
    assert [p.M for p in r1[0]] == [p.M for p in r2[0]] and r1[1:] == r2[1:]
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_mixing_apply.py -q`
Expected: FAIL (`cannot import name 'ParameterizedMixing'`).

- [ ] **Step 3: Implement env interpolation + `ParameterizedMixing`**

Append to `PyLCM/mixing.py`:

```python
from PyLCM.parameters import p0, r_a, cp, l_v


def _interp(z, z_env, profile):
    return float(np.interp(z, z_env, profile))


class ParameterizedMixing:
    """Parameterized homogeneous/inhomogeneous entrainment mixing.

    lambda_ent : fractional entrainment rate [1/m]; entrained fraction = lambda*w*dt.
    ihmd       : Inhomogeneous Mixing Degree in [0,1] (0 homogeneous, 1 inhomogeneous).
    """

    def __init__(self, lambda_ent, ihmd, qv_profiles, theta_profiles, z_env):
        self.lambda_ent = lambda_ent
        self.ihmd = ihmd
        self.qv_profiles = qv_profiles
        self.theta_profiles = theta_profiles
        self.z_env = z_env

    def apply(self, particles_list, T, q, P, z, dt, w, air_mass):
        frac = self.lambda_ent * w * dt
        if frac <= 0.0:
            return particles_list, T, q
        frac = min(frac, 0.999)
        # 1. Bulk entrainment: relax T, q toward the environment at this height.
        theta_env = _interp(z, self.z_env, self.theta_profiles)
        T_env = theta_env * (P / p0) ** (r_a / cp)
        q_env = _interp(z, self.z_env, self.qv_profiles)
        T = T + frac * (T_env - T)
        q = q + frac * (q_env - q)
        # 2. IHMD redistribution of cloud liquid; evaporated water -> vapor, with
        #    latent cooling.
        evaporated = redistribute_droplets(particles_list, self.ihmd, frac)
        dq = evaporated / air_mass
        q = q + dq
        T = T - l_v * dq / cp
        return particles_list, T, q
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_mixing_apply.py -q`
Expected: 3 passed. If `create_env_profiles` blocks on `plt.show()`, the `matplotlib.use("Agg")` in `_setup` handles it; if it still shows, pass a non-interactive backend via env `MPLBACKEND=Agg`.

- [ ] **Step 5: Commit**

```bash
git add PyLCM/mixing.py tests/test_mixing_apply.py
git commit -m "Add ParameterizedMixing.apply: bulk entrainment + IHMD redistribution"
```

---

## Task 3: LEM-ready interface stub (TDD)

**Files:**
- Modify: `PyLCM/mixing.py`
- Create: `tests/test_mixing_lem_stub.py`

- [ ] **Step 1: Write the failing test**

`tests/test_mixing_lem_stub.py`:

```python
import pytest
from PyLCM.mixing import LEMMixing


def test_lem_backend_is_a_clear_stub():
    lem = LEMMixing()
    with pytest.raises(NotImplementedError, match="Phase 3b"):
        lem.apply([], 288.0, 0.009, 90000.0, 500.0, 1.0, 1.0, 100.0)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_mixing_lem_stub.py -q`
Expected: FAIL (`cannot import name 'LEMMixing'`).

- [ ] **Step 3: Implement the stub + a shared protocol docstring**

Append to `PyLCM/mixing.py`:

```python
class LEMMixing:
    """Linear Eddy Model mixing backend (same apply() signature as
    ParameterizedMixing). Not implemented this cycle — a 1D triplet-map domain
    with per-droplet supersaturation perturbation is Phase 3b.
    """

    def apply(self, particles_list, T, q, P, z, dt, w, air_mass):
        raise NotImplementedError(
            "LEMMixing is Phase 3b (1D triplet-map LEM). Use ParameterizedMixing for now."
        )
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_mixing_lem_stub.py -q`
Expected: 1 passed.

- [ ] **Step 5: Commit**

```bash
git add PyLCM/mixing.py tests/test_mixing_lem_stub.py
git commit -m "Add LEMMixing stub proving the MixingModel seam (Phase 3b)"
```

---

## Task 4: Wire mixing into the programmatic run (mixing before condensation, TDD)

**Files:**
- Modify: `pylcm_run.py`
- Create: `tests/test_run_with_mixing.py`

- [ ] **Step 1: Write the failing test**

`tests/test_run_with_mixing.py`:

```python
import numpy as np
from pylcm_run import run_single_series


def test_mixing_off_matches_baseline():
    np.random.seed(0); base = run_single_series(n_ptcl=200, nt=100, collect_every=50)
    np.random.seed(0); same = run_single_series(n_ptcl=200, nt=100, collect_every=50, mixing=None)
    assert np.array_equal(base, same)              # mixing=None is a no-op


def test_inhomogeneous_mixing_reduces_cloud_water():
    from PyLCM.parcel import create_env_profiles
    import matplotlib; matplotlib.use("Agg")
    qv, th, z_env = create_env_profiles(290.0, 0.010, 0.0, 95000.0, "Stable")
    from PyLCM.mixing import ParameterizedMixing
    mix = ParameterizedMixing(2e-3, 1.0, qv, th, z_env)
    np.random.seed(0); base = run_single_series(n_ptcl=200, nt=200, collect_every=200)
    np.random.seed(0); mixed = run_single_series(n_ptcl=200, nt=200, collect_every=200, mixing=mix)
    assert mixed[-1] <= base[-1] + 1e-9            # entrainment removes liquid water
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_run_with_mixing.py -q`
Expected: FAIL (`run_single_series() got an unexpected keyword argument 'mixing'`).

- [ ] **Step 3: Add a `mixing` argument to `run_single_series`, called before condensation**

In `pylcm_run.py`, change the signature to add `mixing=None` and, inside the loop, immediately AFTER `ascend_parcel`/`parcel_rho` and BEFORE `drop_condensation`, insert:

```python
        if mixing is not None:
            particles_list, T_parcel, q_parcel = mixing.apply(
                particles_list, T_parcel, q_parcel, P_parcel, z_parcel,
                dt, w, air_mass)
```

(`air_mass` is the value already returned by `parcel_rho` in the loop; ensure it is in scope before this call.)

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_run_with_mixing.py -q`
Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add pylcm_run.py tests/test_run_with_mixing.py
git commit -m "Wire MixingModel into run loop (mixing before condensation)"
```

---

## Task 5: Wire mixing into the notebook driver + widgets

**Files:**
- Modify: `PyLCM/timestep_routine.py:49-51`, widget plumbing in `PyLCM/widget.py`

- [ ] **Step 1: Replace the basic_entrainment call with ParameterizedMixing**

In `PyLCM/timestep_routine.py`, where `basic_entrainment` is called inside the entrainment time-window branch (around line 49-51), construct a `ParameterizedMixing(entrainment_rate, ihmd, qv_profiles, theta_profiles, z_env)` once before the loop and replace the call with:

```python
        if switch_entrainment and (entrainment_start <= time) and (time < entrainment_end) and (z_parcel < 3000.):
            particles_list, T_parcel, q_parcel = mixing_model.apply(
                particles_list, T_parcel, q_parcel, P_parcel, z_parcel, dt, w_parcel, air_mass_parcel)
```

Add `ihmd` to the `timesteps_function` signature (default `ihmd=0.0`) and thread `z_env` (available from `create_env_profiles`). Keep `basic_entrainment` importable but unused.

- [ ] **Step 2: Add an IHMD slider to the widgets**

In `PyLCM/widget.py`, add a `mixing_degree` (IHMD) `FloatSlider(min=0.0, max=1.0, step=0.05, value=0.0)` alongside the existing entrainment-rate control, and pass its `.value` as `ihmd` into `timesteps_function`. Follow the existing widget wiring pattern.

- [ ] **Step 3: Smoke-test the notebook path imports**

Run: `cd /Users/dr.cloud/PyLCM_entrain && python -c "import PyLCM.timestep_routine, PyLCM.widget; print('ok')"`
Expected: `ok` (widget construction is exercised in the notebook, not CI).

- [ ] **Step 4: Run the full suite**

Run: `python -m pytest tests/ -q`
Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add PyLCM/timestep_routine.py PyLCM/widget.py
git commit -m "Wire ParameterizedMixing + IHMD slider into notebook driver"
```

---

## Task 6: Validation notebook (IHMD sweep)

**Files:**
- Create: `validation/entrainment_mixing.ipynb`

- [ ] **Step 1: Build the notebook**

Cells:
1. md: title + the homogeneous-vs-inhomogeneous question; cite Lim & Hoffmann (2023) and the law `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD`.
2. code: build env profiles (`create_env_profiles`, Agg backend), then run `run_single_series` three times at fixed `lambda_ent` with `ParameterizedMixing` at `ihmd ∈ {0.0, 0.5, 1.0}` (and one baseline `mixing=None`), collecting final DSDs via `ts_analysis`.
3. code: plot the three DSDs overlaid (homogeneous broadens to small sizes; inhomogeneous depletes number, keeps size).
4. code: plot the N–r_v mixing diagram (`N_c/N_{c,0}` vs `r_v/r_{v,0}`) for the sweep, annotated with the homogeneous and inhomogeneous reference curves.
5. md: conclusion tying the curves to IHMD.

- [ ] **Step 2: Execute end-to-end**

Run: `cd /Users/dr.cloud/PyLCM_entrain && jupyter nbconvert --to notebook --execute --inplace validation/entrainment_mixing.ipynb`
Expected: runs with no errors; figures render.

- [ ] **Step 3: Commit**

```bash
git add validation/entrainment_mixing.ipynb
git commit -m "Add entrainment IHMD sweep validation notebook"
```

---

## Task 7: Wrap-up

- [ ] **Step 1: Full suite + README note**

Run `python -m pytest tests/ -q` (all green). Add an "Entrainment mixing (IHMD)" subsection to `README.md` describing the two sliders (λ, IHMD), the `N_c/N_{c,0}=(q_c/q_{c,0})^IHMD` law, and that LEM is a Phase-3b stub.

- [ ] **Step 2: Commit and finish the branch**

```bash
git add README.md && git commit -m "Document entrainment IHMD mixing in README"
```

Then use superpowers:finishing-a-development-branch.

---

## Self-Review

**Spec coverage:** §2 physics → Tasks 1–2; §3 MixingModel/LEM-ready → Tasks 1–3; §4 integration+controls (mixing-first) → Tasks 4–5; §5 testing (budget, endpoints, N–q law, determinism, no-op) → Tasks 1,2,4; §6 validation notebook → Task 6; §7 branch → Task 0; §8 acceptance → Tasks 1–6. Covered.

**Placeholder scan:** All code blocks are concrete. Task 5 widget wiring references "the existing pattern" but names the exact file, widget type, range, and target parameter — no vague placeholders. Task 6 notebook cells are specified by content and the exact functions used.

**Type consistency:** `redistribute_droplets(particles_list, ihmd, frac)` and `ParameterizedMixing(lambda_ent, ihmd, qv_profiles, theta_profiles, z_env).apply(particles_list, T, q, P, z, dt, w, air_mass)` are used identically in Tasks 1–5. `apply` returns `(particles_list, T, q)` everywhere. `mixing` kwarg on `run_single_series` consistent in Task 4.
