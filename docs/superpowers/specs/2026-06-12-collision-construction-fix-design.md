# Integer Multiplicities — Eliminate Zero-Radius Ghost Droplets (SAM-faithful)

**Date:** 2026-06-12
**Author:** J. Lim (jslim93) with agent assistance
**Base:** `feature/entrainment-mixing`
**Status:** Approved approach (Option A — integer multiplicities) — pending spec review → plan

---

## 1. Problem (root-caused empirically)

Collisions can leave a "ghost" super-droplet with **liquid mass `M=0` but multiplicity `A>0`**
(a bare aerosol whose radius the code computes as 0, crashing `E_H80`). Instrumented mechanism
(turbulent maritime, t=971):

1. `Random` init gives every super-droplet an *identical* multiplicity `A`, stored as a **float**.
2. The equal-weight special case routes via **exact** `particle1.A == particle2.A`.
3. After a prior collision two droplets become **near- but not exactly-equal** floats; the `==`
   fails → they fall to the asymmetric `liquid_update_collection`.
4. With `A1 ≈ A2` the loser is collected to `M = M2 − A1·x_int·p_crit = 0.0` (cancellation),
   while `A2 − A1` rounds to `5.8e-11` instead of 0 → a ghost that survives the `A>0` cleanup.

## 2. Why SAM is immune — and the decision

SAM (`micro_coll.f90`) stores the weighting factor `A` as a **`KIND=8` integer**. As a result:
- `A_n .EQ. A_m` is an *exact integer* comparison → **every** equal-weight pair routes to the
  same-weight split (`FLOOR(A/2)`); none ever "almost" matches.
- The asymmetric branches are strict (`A_n .LT. A_m`) and the limiter keeps the loser `A ≥ 1`,
  so the loser is **never fully emptied** → `M = x·A_new > 0` always.
- A `.GT. 1` guard handles the minimum-multiplicity case.

**Decision — Option A: make PyLCM multiplicities integer-valued, matching SAM.** This removes
the entire bug *class* (equality routing, the limiter, and the A≥1 floor all become exact)
rather than mimicking the outcome in float arithmetic.

**Storage choice:** keep `A` in the existing `float64` field but maintain the **invariant that `A`
is always a non-negative integer value**. float64 represents integers exactly up to 2⁵³ ≈ 9e15,
far above the `A ~ 1e10` values here, so integer-valued floats compare and subtract exactly
(`A1 == A2` and `A2 − A1` are exact). This avoids a dtype migration through numba kernels while
giving SAM's exact-integer behavior. (A later cleanup may switch the field to int64.)

## 3. Changes

### 3.0 Parcel-volume scaling (makes the integer conversion lossless)
The parcel volume is physically arbitrary (only concentrations matter), and `A ∝ air_mass`, so
scaling the parcel up lifts every multiplicity — including the rare large-aerosol tail bins —
above 1 without dropping anything. The model is **scale-invariant**: `p_crit = max(A)·K/V·dt`
(both `A` and `V ∝ air_mass`) and the thermodynamics (`ΔM/air_mass`, `ΔM ∝ A ∝ air_mass`) are
unchanged by the scale.
- Replace the two hardcoded `100.0` parcel-air-mass values (`parcel.py parcel_rho`,
  `collision.py determine_collision:191`) with **one shared constant `PARCEL_AIR_MASS`** in
  `parameters.py`, used by both, so `A`'s normalization and collision's `V_parcel` can never drift.
- Set `PARCEL_AIR_MASS` large enough that the smallest kept multiplicity is a comfortable integer
  (target `min A ≳ 100`). A fixed large default (e.g. `1e6` kg) is safe — `A` stays far below the
  `2^53 ≈ 9e15` exact-integer ceiling. (Adaptive sizing in `aero_init` is an alternative if a
  fixed value proves too small/large for extreme aerosol; default to the fixed value.)

### 3.1 Initialization (`aero_init.py`, both modes)
- With the larger `PARCEL_AIR_MASS`, round each `particle.A` to the nearest integer value at
  creation (`A = float(round(A_raw))`) — now **lossless** (all `A` are large integers ≳ 100).
- Drop only super-droplets that still round to `A < 1` (genuinely negligible); with the scaled
  volume this should be empty in normal cases.

### 3.2 Collision (`collision.py`)
- With integer-valued `A`, the existing exact `particle1.A == particle2.A` routing now correctly
  catches all equal-weight pairs → `same_weights_update`.
- `same_weights_update`: split with **floor**, SAM-style: `A_n_new = floor(A/2)`,
  `A_m_new = A − A_n_new`; if either side would be `< 1`, fully merge into one super-droplet
  (the other side's `A → 0`, removed by the end-of-loop `A>0` cleanup) instead of a fractional ghost.
- `liquid_update_collection` (asymmetric): all `A` arithmetic is integer-valued; the existing
  `p_crit` limiter (`A_max − p_crit·A_min ≥ 1`) keeps the loser `A ≥ 1`, so it is never emptied
  and `M2_new = x_int·A2_new > 0`. Compute `M2_new = x_int·A2_new` (reconstruct, not subtract) so
  it is exact. Keep `Ns` consistent the same way.
- Add the SAM `max(A_n, A_m) > 1` guard before coalescence (a single-droplet super-droplet has no
  one to give to).

### 3.3 IHMD entrainment reconciliation (`mixing.py`)
- `redistribute_droplets` currently does `A *= (1−frac)^IHMD` (fractional). Change to **integer
  droplet removal**: `A_new = round(A · (1−frac)^IHMD)` (and `M *= (1−frac)` unchanged). Because
  `PARCEL_AIR_MASS` is now large, every `A` is large, so the discretization error in
  `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD` is negligible.
- Update the redistribution unit tests to use large `A` (e.g. `A=1e8`) so the IHMD law still holds
  to tight tolerance (`rtol≈1e-7`) under integer rounding. Endpoints stay exact (IHMD=0 leaves `A`
  unchanged; IHMD=1 gives `round(A·(1−frac))`).

### 3.4 Defense in depth
Keep the `determine_collision` guard (`M≤0 or A≤0 → skip`) as a backstop; with integer `A` it
should be unreachable for `M=0`.

## 4. Testing

- **No ghosts:** over the turbulent maritime run that failed at t=971, and a 10k-particle ×
  3600-step run, assert **no super-droplet ever has `M≤0` while `A>0`**, no `ZeroDivisionError`,
  and every `A` is integer-valued (`A == round(A)`).
- **Radius positivity:** every `A>0` droplet has `r > 0`.
- **Conservation preserved:** `tests/test_collision_invariants.py`, the condensation golden test,
  and the full suite still pass — total water conservation unchanged.
- **Routing:** exactly-equal integer `A` pairs route to `same_weights_update`; the floor-split
  conserves total `A` and `M`.
- **Both init modes:** ghost-freedom and integer-`A` invariant hold for `Random` and `Weighting_factor`.
- **IHMD law (integer):** `N_c/N_{c,0} ≈ (q_c/q_{c,0})^IHMD` within integer-rounding tolerance.
- **Regression:** re-run `validation/run_physics_checks.py` → still 22/22.

## 5. Acceptance criteria

- Long turbulent + 10k×3600 runs complete with zero `M=0,A>0` droplets, no crash, all `A` integer.
- All existing tests + golden + 22/22 physics checks pass (IHMD law tests adjusted for integer rounding).
- Both init modes produce only physical (r>0), integer-multiplicity droplets.
- No init mode deprecated; default unchanged.

## 6. Branch / workspace

Continue on `feature/entrainment-mixing`. Never touch the live `PyLCM_edu` session. No
"Claude"/Co-Authored-By commit lines.
