# Integer Multiplicities (SAM-faithful) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate zero-radius "ghost" droplets by making super-droplet multiplicity `A` integer-valued (like SAM's `KIND=8` weighting factors), with a large `PARCEL_AIR_MASS` so the integer conversion is lossless.

**Architecture:** `A` stays in the `float64` field but always holds an exact integer. A single `PARCEL_AIR_MASS` constant (large) replaces the hardcoded `100.0` parcel mass in `parcel_rho` and `determine_collision`, scaling all `A` to large integers. Collision uses integer-exact equal-weight routing + floor-split + reconstructed `M=x·A`; IHMD mixing removes whole droplets. The model is scale-invariant, so intensive quantities (T, q, concentrations) are unchanged.

**Tech Stack:** Python, numpy, pytest. Worktree `/Users/dr.cloud/PyLCM_entrain` (branch `feature/entrainment-mixing`).

**Constraints:** Never add a "Claude"/Co-Authored-By commit line. Never touch `/Users/dr.cloud/PyLCM_edu`. Run with `MPLBACKEND=Agg`.

---

## File Structure

| Path | Change |
|------|--------|
| `PyLCM/parameters.py` | Add `PARCEL_AIR_MASS` constant |
| `PyLCM/parcel.py` | `parcel_rho` uses `PARCEL_AIR_MASS` instead of `100.0` |
| `PyLCM/collision.py` | `determine_collision` uses `PARCEL_AIR_MASS`; integer construction in `liquid_update_collection`, `same_weights_update`; `max(A)>1` guard |
| `PyLCM/aero_init.py` | Round `A` to integer in both init modes; drop `A<1` |
| `PyLCM/mixing.py` | `redistribute_droplets` integer droplet removal |
| `tests/` | New ghost-freedom + integer-A tests; updated IHMD-law + golden tests |
| `validation/golden_setup.py`, `tests/golden/condensation_golden.npz` | Regenerate golden under the new scale |

---

## Task 0: Baseline

- [ ] **Step 1: Confirm starting state**

Run: `cd /Users/dr.cloud/PyLCM_entrain && MPLBACKEND=Agg python -m pytest tests/ -q`
Expected: `33 passed`. We are on branch `feature/entrainment-mixing`.

---

## Task 1: `PARCEL_AIR_MASS` constant, unified and scale-invariant

**Files:** Modify `PyLCM/parameters.py`, `PyLCM/parcel.py`, `PyLCM/collision.py`; Create `tests/test_scale_invariance.py`

- [ ] **Step 1: Write the scale-invariance test (fails until the constant exists & is used)**

`tests/test_scale_invariance.py`:

```python
import numpy as np
import matplotlib; matplotlib.use("Agg")
from PyLCM import parameters
from PyLCM.parcel import parcel_rho


def test_air_mass_constant_drives_parcel_rho():
    rho, V, air_mass = parcel_rho(95000.0, 285.0)
    assert abs(air_mass - parameters.PARCEL_AIR_MASS) / parameters.PARCEL_AIR_MASS < 1e-9
    assert abs(V - parameters.PARCEL_AIR_MASS / rho) / V < 1e-9
```

- [ ] **Step 2: Run to verify it fails**

Run: `MPLBACKEND=Agg python -m pytest tests/test_scale_invariance.py -q`
Expected: FAIL (`module 'PyLCM.parameters' has no attribute 'PARCEL_AIR_MASS'`).

- [ ] **Step 3: Add the constant** to `PyLCM/parameters.py` (after the `pi` line):

```python
# Parcel air mass (kg). Physically arbitrary in a parcel model (only concentrations
# matter); chosen large so super-droplet multiplicities A = air_mass*concentration are
# large integers (lossless integer rounding, no dropped bins). Used by parcel_rho and
# the collision volume so the two can never drift.
PARCEL_AIR_MASS = 1.0e6
```

- [ ] **Step 4: Use it in `parcel_rho`** (`PyLCM/parcel.py`) — replace the `100.0` line:

```python
    V_parcel   = PARCEL_AIR_MASS / rho_parcel # volume of the parcel for PARCEL_AIR_MASS of air
    air_mass_parcel = V_parcel * rho_parcel
```

Ensure `PyLCM/parcel.py` imports it: it already does `from PyLCM.parameters import *`, so `PARCEL_AIR_MASS` is in scope.

- [ ] **Step 5: Use it in `determine_collision`** (`PyLCM/collision.py:190-191`) — replace `V_parcel = 100.0 / rho_parcel`:

```python
    # V_parcel = air_mass / rho_parcel, consistent with parcel_rho (PARCEL_AIR_MASS)
    V_parcel = PARCEL_AIR_MASS / rho_parcel
```

`collision.py` imports parameters via `from PyLCM.parameters import *` (verify; if not, add `from PyLCM.parameters import PARCEL_AIR_MASS`).

- [ ] **Step 6: Run the scale-invariance + full suite**

Run: `MPLBACKEND=Agg python -m pytest tests/test_scale_invariance.py tests/test_collision_invariants.py -q`
Expected: scale test PASS. The collision invariant tests still pass (they construct particles directly, unaffected). The golden condensation test will now FAIL (absolute M/A scaled) — that is expected and fixed in Task 5; do not "fix" it here.

- [ ] **Step 7: Commit**

```bash
git add PyLCM/parameters.py PyLCM/parcel.py PyLCM/collision.py tests/test_scale_invariance.py
git commit -m "Add unified PARCEL_AIR_MASS constant (large parcel volume) for parcel_rho and collision"
```

---

## Task 2: Integer multiplicities at initialization

**Files:** Modify `PyLCM/aero_init.py:178` (Random) and `:208` (Weighting_factor); Create `tests/test_integer_multiplicity.py`

- [ ] **Step 1: Write the failing test**

`tests/test_integer_multiplicity.py`:

```python
import numpy as np
from PyLCM.aero_init import aero_init


def _init(mode):
    mu = np.log(np.array([0.02e-6, 0.2e-6])); sigma = np.log(np.array([1.4, 1.6]))
    return aero_init(mode, 200, 95000.0, 0.0, 285.0, 0.008,
                     np.array([100e6, 20e6]), mu, sigma, 1777.0,
                     np.array([0.5, 0.5]), True)[2]


def test_random_init_multiplicities_are_integers_ge_1():
    pl = _init("Random")
    assert pl and all(p.A == round(p.A) and p.A >= 1 for p in pl)


def test_weighting_factor_init_multiplicities_are_integers_ge_1():
    pl = _init("Weighting_factor")
    assert pl and all(p.A == round(p.A) and p.A >= 1 for p in pl)
```

- [ ] **Step 2: Run to verify it fails**

Run: `MPLBACKEND=Agg python -m pytest tests/test_integer_multiplicity.py -q`
Expected: FAIL (float `A` not integer-valued).

- [ ] **Step 3: Round `A` in Random mode** (`aero_init.py:178`) — replace `particle.A = air_mass_parcel * np.sum(N_aero)/n_ptcl`:

```python
            particle.A = float(round(air_mass_parcel * np.sum(N_aero) / n_ptcl))
```

- [ ] **Step 4: Round `A` in Weighting_factor mode** (`aero_init.py:208`) — replace `particle.A = air_mass_parcel * pdf_sum[i] * dlogr * radius[i]`:

```python
            particle.A = float(round(air_mass_parcel * pdf_sum[i] * dlogr * radius[i]))
```

- [ ] **Step 5: Drop sub-1 (A<1) super-droplets after the init loop** — at the end of `aero_init`, before the `dql_liq` recomputation (line ~224), insert:

```python
    # Remove super-droplets whose rounded multiplicity is below one real droplet.
    particles_list = [p for p in particles_list if p.A >= 1]
```

- [ ] **Step 6: Run to verify it passes**

Run: `MPLBACKEND=Agg python -m pytest tests/test_integer_multiplicity.py -q`
Expected: 2 passed.

- [ ] **Step 7: Commit**

```bash
git add PyLCM/aero_init.py tests/test_integer_multiplicity.py
git commit -m "Round multiplicities to integers at init (both modes); drop A<1"
```

---

## Task 3: Integer-exact collision construction (no ghosts)

**Files:** Modify `PyLCM/collision.py` (`collection`, `liquid_update_collection`, `same_weights_update`); Create `tests/test_no_ghost_droplets.py`

- [ ] **Step 1: Write the ghost-freedom test**

`tests/test_no_ghost_droplets.py`:

```python
import numpy as np
import matplotlib; matplotlib.use("Agg")
from validation.phys_harness import run


def test_turbulent_collision_makes_no_ghosts():
    # The config that previously crashed at t~971 with a zero-radius ghost.
    _, pl = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1800,
                collisions=True, switch_turb=True, eps=0.04)
    # No super-droplet may have zero/negative liquid mass while keeping A>0,
    # and every surviving A must be a positive integer.
    for p in pl:
        assert p.A == round(p.A)
        if p.A > 0:
            assert p.M > 0.0
```

- [ ] **Step 2: Run to verify it fails (or errors)**

Run: `MPLBACKEND=Agg python -m pytest tests/test_no_ghost_droplets.py -q`
Expected: FAIL — a droplet with `M==0` and `A>0` exists (ghost), or `A` non-integer.

- [ ] **Step 3: Reconstruct `M` (and `Ns`) in `liquid_update_collection`** (`collision.py:106-109`) — replace the two subtraction lines for the loser with reconstruction from the updated multiplicity:

```python
    #Decrease of number, aerosol and water mass due to collision
    ptcl_int2.A  = ptcl_int2.A - ptcl_int1.A * p_crit
    # Reconstruct from per-droplet mass so M, Ns stay exactly proportional to A
    # (avoids catastrophic cancellation that left M=0 while A>0).
    ptcl_int2.M  = x_int * ptcl_int2.A
    ptcl_int2.Ns = xs_int * ptcl_int2.A
```

- [ ] **Step 4: Floor-split in `same_weights_update`** (`collision.py:159-165`) — replace the `A*0.5` split with an integer floor split, merging fully if a half would be < 1:

```python
    # SAM-style integer floor split of equal weighting factors.
    A_total = ptcl_int1.A
    A_half = float(A_total // 2)
    if A_half < 1:
        # Too few to split: merge all into ptcl_int1, remove ptcl_int2.
        ptcl_int1.A  = A_total
        ptcl_int1.M  = (xn + xm) * A_total
        ptcl_int1.Ns = (xsn + xsm) * A_total
        ptcl_int2.A  = 0.0
        ptcl_int2.M  = 0.0
        ptcl_int2.Ns = 0.0
    else:
        ptcl_int1.A  = A_half
        ptcl_int2.A  = A_total - A_half
        ptcl_int1.M  = (xn + xm) * ptcl_int1.A
        ptcl_int2.M  = (xn + xm) * ptcl_int2.A
        ptcl_int1.Ns = (xsn + xsm) * ptcl_int1.A
        ptcl_int2.Ns = (xsn + xsm) * ptcl_int2.A
```

(Note: `same_weights_update` is only called when `particle1.A == particle2.A`, so `ptcl_int1.A == ptcl_int2.A == A_total` on entry; the existing `xn,xm,xsn,xsm` are computed above this block — keep them.)

- [ ] **Step 5: Add the `max(A) > 1` coalescence guard** in `collection()` — after the existing `if min(particle1.A, particle2.A) <= 0: continue` (line 35), add:

```python
        # A single-droplet super-droplet (A<=1) has no droplets to give away.
        if max(particle1.A, particle2.A) <= 1:
            continue
```

- [ ] **Step 6: Run ghost-freedom + collision invariants**

Run: `MPLBACKEND=Agg python -m pytest tests/test_no_ghost_droplets.py tests/test_collision_invariants.py tests/test_collision_zero_radius.py -q`
Expected: all pass. (`test_collision_invariants` proves mass conservation still holds with reconstruction.)

- [ ] **Step 7: Commit**

```bash
git add PyLCM/collision.py tests/test_no_ghost_droplets.py
git commit -m "Integer-exact collision: reconstruct M=x*A, floor-split equal weights, A>1 guard"
```

---

## Task 4: Integer droplet removal in IHMD mixing

**Files:** Modify `PyLCM/mixing.py` (`redistribute_droplets`); Modify `tests/test_mixing_redistribute.py`

- [ ] **Step 1: Update the redistribution tests for integer removal (large A)**

In `tests/test_mixing_redistribute.py`, change `_cloud` to use a large integer `A` and relax the law tolerance:

```python
def _cloud(n=50, M=2.0e-9, A=1.0e8):
    out = []
    for _ in range(n):
        p = particles(1); p.M, p.A, p.Ns, p.kappa = M, A, 1e-18, 0.5
        out.append(p)
    return out
```

In `test_nq_power_law_holds_exactly` change the tolerance and name:

```python
def test_nq_power_law_holds_under_integer_removal():
    for ihmd in (0.0, 0.5, 1.0):
        ps = _cloud()
        N0 = sum(p.A for p in ps); q0 = sum(p.M for p in ps)
        evap = redistribute_droplets(ps, ihmd=ihmd, frac=0.3)
        N1 = sum(p.A for p in ps); q1 = sum(p.M for p in ps)
        assert all(p.A == round(p.A) for p in ps)             # integer multiplicities
        assert abs(evap - 0.3 * q0) < 1e-18                   # evaporated mass exact
        assert np.isclose(N1 / N0, (q1 / q0) ** ihmd, rtol=1e-6)  # IHMD law (integer)
```

In `test_homogeneous_conserves_number`, the per-droplet mass assertion uses `2.0e-9/100.0`; update the constant to match `A=1e8`: `(2.0e-9/1.0e8)*0.6`. In `test_inhomogeneous_preserves_droplet_size`, update `2.0e-9/100.0` to `2.0e-9/1.0e8`.

- [ ] **Step 2: Run to verify it fails**

Run: `MPLBACKEND=Agg python -m pytest tests/test_mixing_redistribute.py -q`
Expected: FAIL (current `A *= keep_num` yields non-integer `A`).

- [ ] **Step 3: Integer removal in `redistribute_droplets`** (`mixing.py`) — replace `p.A = p.A * keep_num`:

```python
        p.M = p.M * keep_mass
        p.A = float(round(p.A * keep_num))
```

(Keep the `M *= keep_mass` as the liquid reduction; only `A` becomes integer.)

- [ ] **Step 4: Run to verify it passes**

Run: `MPLBACKEND=Agg python -m pytest tests/test_mixing_redistribute.py tests/test_mixing_apply.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add PyLCM/mixing.py tests/test_mixing_redistribute.py
git commit -m "IHMD redistribution removes whole droplets (integer A)"
```

---

## Task 5: Regenerate golden; verify scale-invariance; full validation

**Files:** Modify `tests/golden/condensation_golden.npz` (regenerate); Create `tests/test_scale_invariant_intensive.py`

- [ ] **Step 1: Capture the OLD golden's intensive quantities, then prove invariance**

`tests/test_scale_invariant_intensive.py`:

```python
import numpy as np
from validation.golden_setup import run_condensation_only

# Recorded from the pre-change run (air_mass=100, float A): condensation is
# scale-invariant, so T and q must be unchanged (to integer-rounding tolerance)
# after the PARCEL_AIR_MASS + integer-A change.
OLD_T = 291.302817
OLD_Q = 0.01358930


def test_condensation_T_q_are_scale_invariant():
    out = run_condensation_only()
    assert np.isclose(out["T"], OLD_T, rtol=1e-4)
    assert np.isclose(out["q"], OLD_Q, rtol=1e-3)
```

- [ ] **Step 2: Run it — confirms physics unchanged by the scaling**

Run: `MPLBACKEND=Agg python -m pytest tests/test_scale_invariant_intensive.py -q`
Expected: PASS (T, q within tolerance of the pre-change values — the change is scale-invariant).
If this FAILS, the change altered intensive physics — STOP and investigate before regenerating the golden.

- [ ] **Step 3: Regenerate the golden snapshot under the new scale**

Run: `MPLBACKEND=Agg python -m validation.make_golden`
Expected: prints `Golden saved: …`. This overwrites `tests/golden/condensation_golden.npz` with the new (scaled, integer-A) absolute M/A and the (unchanged) T/q.

- [ ] **Step 4: Confirm the golden test passes against the regenerated snapshot**

Run: `MPLBACKEND=Agg python -m pytest tests/test_golden_condensation.py tests/test_condensation_fast.py -q`
Expected: 4 passed (the fast SoA path must still bit-match the regenerated golden).

- [ ] **Step 5: Full suite + cross-feature physics validation**

Run: `MPLBACKEND=Agg python -m pytest tests/ -q`
Expected: all pass.
Run: `MPLBACKEND=Agg python -m validation.run_physics_checks 2>&1 | tail -3`
Expected: `22/22 physical checks PASS`.

- [ ] **Step 6: Commit**

```bash
git add tests/golden/condensation_golden.npz tests/test_scale_invariant_intensive.py
git commit -m "Regenerate condensation golden under PARCEL_AIR_MASS scale; assert T/q scale-invariance"
```

---

## Task 6: Robustness re-check at scale and wrap-up

- [ ] **Step 1: 10k x 3600 ghost-free + no-crash check (the original ensemble crash)**

Run:
```bash
cd /Users/dr.cloud/PyLCM_entrain && MPLBACKEND=Agg python -c "
import warnings; warnings.filterwarnings('ignore')
from validation.phys_harness import run
_, pl = run(seed=0, aerosol='maritime', n_ptcl=10000, nt=3600, collisions=True, switch_turb=True, eps=0.04)
bad = [p for p in pl if p.A>0 and p.M<=0]
assert not bad and all(p.A==round(p.A) for p in pl), bad
print('10k x 3600 turbulent: no ghosts, all integer A, no crash')
"
```
Expected: prints the success line (this is the configuration that previously raised `ZeroDivisionError`).

- [ ] **Step 2: Update README note**

In `README.md`, under the collision/validation text, add one line: multiplicities are integer-valued (SAM-faithful) with a large `PARCEL_AIR_MASS`, eliminating zero-radius ghost droplets.

- [ ] **Step 3: Commit and finish**

```bash
git add README.md
git commit -m "Document integer-multiplicity collision fix"
```

Then use superpowers:finishing-a-development-branch.

---

## Self-Review

**Spec coverage:** §3.0 PARCEL_AIR_MASS→Task 1; §3.1 integer init→Task 2; §3.2 collision construction→Task 3; §3.3 IHMD integer→Task 4; §3.4 backstop guard (kept, unchanged); §4 testing→Tasks 1-6 (ghost-freedom, integer-A, IHMD-law, golden, 22/22). Covered.

**Placeholder scan:** All steps contain concrete code/commands. Golden regeneration (Task 5) intentionally overwrites a binary artifact via the existing `make_golden` script — not a placeholder.

**Type consistency:** `PARCEL_AIR_MASS` used identically in `parameters.py`/`parcel.py`/`collision.py`. `A = float(round(...))` pattern consistent across `aero_init` (Task 2), `same_weights_update`/`liquid_update_collection` (Task 3), `redistribute_droplets` (Task 4). `run(...)` harness signature matches `validation/phys_harness.py`. Golden keys `M/A/T/q` consistent with `golden_setup.py`.
