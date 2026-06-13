# Collision Construction Fix — Eliminate Zero-Radius Ghost Droplets

**Date:** 2026-06-12
**Author:** J. Lim (jslim93) with agent assistance
**Base:** `feature/entrainment-mixing` (continues the zero-radius work already there)
**Status:** Approved approach (Path A) — pending spec review → implementation plan

---

## 1. Problem (root-caused empirically)

Collisions can leave a "ghost" super-droplet with **liquid mass `M=0` but multiplicity `A>0`
and aerosol `Ns>0`** — physically a bare aerosol with a real haze radius (~0.29 µm), but the
code computes its radius from liquid mass only → `r=0`, which crashes `E_H80`
(`min(r1/r2, r2/r1)` divides by zero) and is unphysical.

**Exact mechanism (instrumented, turbulent maritime case, t=971):**
1. `Random` init gives every super-droplet identical `A` (Shima 2009 equal-multiplicity, the
   SAM/Fortran-reference setup), so equal-weight collisions are common.
2. `collection()` routes equal weights to `same_weights_update` via **exact** `particle1.A == particle2.A`.
3. After a prior collision, two droplets become **near- but not exactly-equal** in `A`; the `==`
   check fails → they fall to the asymmetric `liquid_update_collection`.
4. With `A1 ≈ A2`, the loser transfers essentially all water: computed
   `M2_new = M2 − A1·x_int·p_crit = M2 − M2 = 0.0` (catastrophic cancellation), while
   `A2_new = A2 − A1` rounds to a tiny positive (`5.8e-11`) instead of 0.
5. End-of-loop cleanup removes only `A ≤ 0`, so the `M=0, A=5.8e-11` ghost **survives** and
   crashes the next collision it enters.

Note: although `Random` triggers it most, near-equal weights can also arise mid-run in
`Weighting_factor` mode, so the fix is in the **shared construction**, not in one init mode.

## 2. Decision

**Path A: keep both init modes; fix the shared `collection` construction.** Do not deprecate
`Random` (it is the equal-multiplicity SDM standard and matches the Fortran reference). The
default init mode is unchanged this cycle (everything downstream uses `Random`); switching the
educational default to `Weighting_factor` for determinism is a separate, later choice.

## 3. The fix (three linked parts, all in `PyLCM/collision.py`)

1. **Tolerant equal-weight routing.** In `collection()`, replace the exact
   `particle1.A == particle2.A` with a relative-tolerance check (e.g. `np.isclose(A1, A2,
   rtol=1e-9)`), so near-equal pairs use `same_weights_update` (the proper equal-weight split)
   instead of the asymmetric path that fully collects the loser.

2. **Reconstruct, don't subtract; snap a fully-collected droplet to zero.** In
   `liquid_update_collection`, compute the loser's new mass as `M2_new = x_int · A2_new`
   (algebraically identical, no catastrophic cancellation), and when `A2_new` is essentially
   fully collected (within a small relative epsilon of 0), set `A2 = 0` so the existing
   end-of-loop cleanup removes it cleanly. No `M=0, A>0` ghost can persist.

3. **Wet-aerosol floor on `M`.** In `liquid_update_collection` and `same_weights_update`, floor
   every surviving droplet's `M` at its wet-aerosol-equilibrium mass — the same floor
   `aero_init` uses:
   `M_floor = max(r_aero, r_equi(S, T, r_aero, κ))**3 · A · 4/3·π·ρ_liq`,
   with `r_aero = (Ns/(A·4/3·π·ρ_aero))**(1/3)`. This guarantees `r ≥ r_equi > 0` by
   construction. Requires threading the parcel supersaturation `S` (or `q`,`T`,`P` to compute
   it) and `T` into `collection()` → the update functions; import `r_equi` locally to avoid a
   circular import.

## 4. `same_weights_update` hardening

- Apply the same **wet-aerosol floor** (§3.3) to both outputs.
- Guard the **small/odd-multiplicity** special case: the current `A*0.5` split can produce a
  sub-1 (fractional) multiplicity. Keep the split's total `A` and `M` conserved, but ensure no
  output has `0 < A < ` a physical minimum that would later misbehave; for `A` below the
  minimum, fully merge into one super-droplet (the other branch's `A → 0`, removed by cleanup)
  rather than emitting a fractional ghost. (Exact threshold finalized in the plan.)

## 5. Defense in depth

Keep the existing `determine_collision` guard (`M≤0 or A≤0 → no collision`) as a cheap backstop,
but it should now be effectively unreachable for `M=0` because the construction no longer
produces such states.

## 6. Testing

- **No ghosts:** over a long turbulent maritime run (the case that failed at t=971) and a
  10k-particle / 3600-step run, assert that **no super-droplet ever has `M≤0` while `A>0`**, and
  no `ZeroDivisionError`.
- **Radius positivity:** every super-droplet with `A>0` has `r ≥ r_equi > 0` at all times.
- **Conservation preserved:** the existing collision invariants (`tests/test_collision_invariants.py`),
  the condensation golden test, and the full suite still pass — the construction change must not
  alter total water conservation.
- **Routing:** two near-equal-`A` droplets route to `same_weights_update`, not the asymmetric path.
- **Both modes:** ghost-freedom holds for `Random` and `Weighting_factor` init.
- **Regression:** re-run `validation/run_physics_checks.py` → still 22/22.

## 7. Acceptance criteria

- Long turbulent + 10k×3600 runs complete with zero `M=0,A>0` droplets and no crash.
- All existing tests + the golden condensation test pass; 22/22 physics checks pass.
- `Random` and `Weighting_factor` both produce only physical (r>0) droplets.
- No init mode deprecated; default unchanged.

## 8. Branch / workspace

Continue on `feature/entrainment-mixing` (this directly extends the zero-radius collision work
already committed there), or a child branch off it. Never touch the live `PyLCM_edu` session.
No "Claude"/Co-Authored-By commit lines.
