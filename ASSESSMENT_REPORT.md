# PyLCM_edu Codebase Assessment Report

**Date:** 2026-02-09
**Branch:** `feature/lecture-readiness-improvements`
**Purpose:** Comprehensive assessment for graduate student lecture readiness
**Priority Order:** 1) Correct Science, 2) Educational Benefit, 3) Code Speed/Quality

---

## Executive Summary

PyLCM_edu is a well-designed Lagrangian Cloud Model educational simulator with solid warm-rain physics (condensation, collision-coalescence, Beard terminal velocity, Hall collision efficiencies). However, the assessment team identified **6 critical science bugs**, **significant pedagogical gaps**, and **several performance issues** that must be addressed before classroom deployment.

**Overall Readiness Scores:**
| Dimension | Score | Status |
|-----------|-------|--------|
| Science Correctness | 6/10 | Collision module has fundamental bugs; condensation is correct |
| Educational Value | 5/10 | Physics engine is strong, but pedagogical wrapper is absent |
| Code Quality | 4/10 | No tests, no error handling, wildcard imports |
| UI/UX | 6/10 | Functional but memory leak in animation |

---

## PRIORITY 1: CORRECT SCIENCE (Must-Fix Before Any Lecture)

### ~~CRITICAL-1: Surface Tension Formula~~ — VERIFIED CORRECT
**File:** `condensation.py:97-103`
**Status:** Initially flagged as wrong, but **verified against the reference Fortran LCM code** (`SAM6.10.10.LCM_JS/SRC/MICRO_LAGRANGE/micro_cond.f90:1475-1481`). The polynomial with all positive coefficients **IS the correct P&K (1997) Eq. 5-12 form**. Both Fortran and Python implementations are identical. No fix needed.

---

### CRITICAL-1 (renumbered): Collision Probability Off by ~82x (V_parcel and rho_parcel Hardcoded)
**File:** `collision.py:173-179`
**Impact:** ALL collision/coalescence results are quantitatively wrong

The `determine_collision()` function receives `rho_parcel` and could compute `V_parcel`, but immediately overwrites them:
```python
rho_parcel = 1.0    # line 175: OVERWRITES actual density (~1.2 kg/m3)
V_parcel = 1.0      # line 178: OVERWRITES actual volume (~82 m3)
```

The weighting factors `A` are computed for a 100 kg parcel (V ~ 82 m3), but collision probability divides by V_parcel = 1.0, making collisions ~82x too frequent. Terminal velocities also use wrong air density.

**Fix:** Remove hardcoded values. Use the function parameter `rho_parcel` and compute `V_parcel = 100.0 / rho_parcel`.

---

### CRITICAL-3: Hall (1980) Collision Efficiency Interpolation Uses Wrong Units
**File:** `collision.py:269`
**Impact:** Collision efficiency interpolation always returns lower-bound values

The interpolation formula uses `rmax` in meters, but the lookup table `r0` is in microns:
```python
pp = (rmax - r0[ir-1]) / (r0[ir] - r0[ir-1])  # rmax in meters, r0 in microns!
```

`rmax ~ 50e-6 m` vs `r0 ~ 50.0 um` makes `pp` essentially zero, disabling interpolation.

**Fix:** Use `rmax * 1.0E6` in the interpolation.

---

### CRITICAL-4: Kappa Corruption in Same-Weight Collision
**File:** `collision.py:162-164`
**Impact:** Hygroscopicity parameter silently corrupted after equal-weight collisions

```python
ptcl_int1.kappa = (v1*ptcl_int1.kappa + v2*ptcl_int2.kappa) / (v1+v2)  # modifies ptcl_int1
ptcl_int2.kappa = (v1*ptcl_int1.kappa + v2*ptcl_int2.kappa) / (v1+v2)  # uses MODIFIED ptcl_int1!
```

Line 163 changes `ptcl_int1.kappa` before line 164 reads it. `ptcl_int2` gets a double-weighted value.

**Fix:** `new_kappa = ...; ptcl_int1.kappa = new_kappa; ptcl_int2.kappa = new_kappa`

---

### CRITICAL-5: N_l Is a Tuple, Not a Float
**File:** `parameters.py:28`

```python
N_l = N_avogadro * rho_liq / molecular_weight_water,  # trailing comma!
```

The trailing comma makes `N_l` a single-element tuple `(3.34e28,)`. Not currently used but is a landmine.

---

### CRITICAL-6: Division by Zero in DSD Analysis
**File:** `analysis.py:67-68`
**Impact:** Crash when no cloud droplets exist (initialization, early timesteps)

```python
rc_liq_avg = np.nansum(...) / np.sum(np.array(particles_ac))  # zero denominator!
```

When all particles are sub-activation-radius, the denominator is zero.

---

### IMPORTANT Science Issues

| Issue | File | Description |
|-------|------|-------------|
| Temperature sign logic | `parcel.py:39-42` | Checks `w_parcel` sign instead of `dz` sign. Correct by accident for w>0, wrong for w<0. Should use `T -= dz*g/cp` unconditionally |
| Inconsistent e_sat | `aero_init.py:124` | Uses Magnus formula; `condensation.py` uses Flatau polynomial. Small (~0.5%) discrepancy at first timestep |
| Wrong V_parcel formula | `condensation.py:15` | `100/P/(R*T)` instead of `100*R*T/P`. Dead code (unused) but misleading |
| No ventilation | `condensation.py:31` | `f_vent = 1.0` hardcoded. Acceptable simplification but should be documented |
| No radiation | `condensation.py:33` | `D_pre = 0.0` hardcoded. Acceptable for educational purposes |
| Weighting_factor mode ignores per-mode kappa | `aero_init.py:211` | All particles get `kappa=0.5` regardless of mode |
| Entrainment doesn't dilute particles | `entrainment.py` | Only modifies T and q, does not reduce droplet number/mass |
| rho_a in wrong units | `parameters.py:3` | g/cm3 while everything else is SI kg/m3 |

---

## PRIORITY 2: EDUCATIONAL BENEFIT (High-Impact Improvements)

### Current State: 5/10 as a teaching tool (8/10 as a research demo)

The notebook is a *control panel* (widgets -> run -> plots), not a *learning experience*. There are only 6 markdown cells, all headers. Zero explanatory text, zero equations, zero guided experiments.

### Week 1: Essential for Lecture Use

1. **Add conceptual markdown cells** throughout the notebook
   - What is a Lagrangian cloud model? Why super-droplets?
   - Kohler theory with LaTeX equations: `S = a/r - b*r_N^3/r^3`
   - Diffusional growth equation with physical explanation
   - Collection kernel and LSM explanation
   - Learning objectives at the top

2. **Create 4 guided experiments** (separate notebook sections):
   - **Experiment 1 — Condensation Only:** Compare supersaturation for w=0.5, 1.0, 2.0 m/s
   - **Experiment 2 — Maritime vs Continental:** Change aerosol modes, compare activation
   - **Experiment 3 — Collision Effects:** Enable collision, observe DSD broadening
   - **Experiment 4 — Full Microphysics:** Design a realistic cumulus simulation

3. **Add standalone Kohler curve visualization** cell before the simulation
   - Plot equilibrium curves for each aerosol mode
   - Show critical supersaturation and critical radius

4. **Add learning objectives** at the top of the notebook

### Week 2: Significantly Enhances Teaching

5. **Add parameter presets** (Maritime/Continental/Arctic auto-fill buttons)
6. **Add annotations to plots** (cloud base height, collision onset, max S markers)
7. **Add reference lines** (dry/moist adiabatic lapse rates on T-z plot)
8. **Create Quick Reference cell** with variable definitions and units
9. **Add self-check questions** after each experiment

### Weeks 3-4: Polish and Depth

10. **Add run comparison capability** (overlay two simulations)
11. **Create homework problem sets** as separate cells
12. **Add references** to textbooks (Rogers & Yau, Pruppacher & Klett)
13. **Create instructor notes** with expected results and common misconceptions

---

## PRIORITY 3: CODE SPEED & QUALITY

### Performance (most impactful first)

| Fix | File | Impact | Effort |
|-----|------|--------|--------|
| Replace `fig.show()` with `FigureWidget` + `batch_update()` | `animation.py` | Eliminates memory leak, 10x rendering speedup | 2-3 hrs |
| Add `@jit(nopython=True)` to `radius_liquid_euler()` | `condensation.py` | 10-50x condensation speedup | 1-2 hrs |
| Vectorize `ts_analysis()` | `analysis.py` | Faster DSD computation | 2-3 hrs |
| Fix `deepcopy` inside loop | `print_plot.py:159` | Faster plot rendering | 10 min |
| Move `time_array` creation outside loop | `timestep_routine.py:81` | Minor speedup | 5 min |

### Code Quality

| Fix | File | Impact | Effort |
|-----|------|--------|--------|
| Change `nt_widget` to `BoundedIntText(min=1)` | `widget.py` | Prevent kernel freeze | 5 min |
| Set sigma min=0.01, mu min=0.001 in grid | `widget.py` | Prevent NaN/div-by-zero | 5 min |
| Add `max_iter=10000` guard to `r_equi()` | `aero_init.py:266` | Prevent infinite loop | 10 min |
| Replace 'jet' colormap with 'viridis' | `print_plot.py` | Scientific best practice | 5 min |
| Fix "Evaporaton" typo | `print_plot.py:194,203` | Professionalism | 1 min |
| Add `tqdm` progress bar to timestep loop | `timestep_routine.py` | Major UX improvement | 30 min |
| Replace `0.33333333333` with `(1.0/3.0)` | Multiple files | Full float64 precision | 15 min |
| Remove unused imports (numba, itertools, lognorm, pylab) | Multiple files | Clean dependencies | 10 min |
| Replace deprecated `cm.get_cmap('jet')` | `print_plot.py` | Compatibility | 5 min |
| Add "Stop" button for simulation | `widget.py` | UX safety | 1-2 hrs |
| Replace wildcard `import *` with explicit imports | All files | Maintainability | 3-4 hrs |

---

## Implementation Roadmap (4-Week Plan)

### Week 1: Critical Science Fixes + Minimal Pedagogy
- [x] ~~Fix surface tension formula~~ — Verified correct against Fortran reference
- [ ] Fix V_parcel/rho_parcel hardcoded in collision (CRITICAL-1)
- [ ] Fix E_H80 interpolation units — missing *1.0E6 conversion (CRITICAL-2)
- [ ] Fix kappa corruption in same_weights_update (CRITICAL-3)
- [ ] Fix N_l trailing comma tuple (CRITICAL-4)
- [ ] Fix division by zero in analysis (CRITICAL-5)
- [ ] Fix temperature sign logic in parcel.py
- [ ] Add learning objectives and intro markdown to notebook
- [ ] Add basic Kohler curve visualization cell

### Week 2: Educational Content + UI Fixes
- [ ] Write 4 guided experiments with instructions
- [ ] Add parameter presets (Maritime/Continental/Arctic)
- [ ] Add conceptual explanations before each widget section
- [ ] Replace Plotly fig.show() with FigureWidget
- [ ] Add tqdm progress bar
- [ ] Fix nt_widget bounds, sigma/mu bounds

### Week 3: Visualization + Performance
- [ ] Add plot annotations (cloud base, collision onset)
- [ ] Replace jet colormap with viridis
- [ ] Add @jit to radius_liquid_euler
- [ ] Vectorize ts_analysis
- [ ] Add self-check questions after experiments

### Week 4: Polish + Testing
- [ ] Add homework problem sets
- [ ] Add Quick Reference cell
- [ ] Create instructor notes
- [ ] Remove unused imports, replace wildcard imports
- [ ] Manual testing of all guided experiments
- [ ] Fix remaining minor issues

---

## Fortran Reference Verification

Compared against the parent Fortran LCM code (`SAM6.10.10.LCM_JS/SRC/MICRO_LAGRANGE/`):

| Item | Python vs Fortran | Status |
|------|------------------|--------|
| Surface tension (P&K Eq 5-12) | **Identical** coefficients and signs | CORRECT |
| Newton-Raphson condensation | **Identical** structure and logic | CORRECT |
| V_parcel in collision probability | Fortran uses `dV(k)` (actual volume), Python uses `1.0` | **BUG CONFIRMED** |
| E_H80 radius interpolation | Fortran converts `r*1.0E6` to microns, Python uses meters | **BUG CONFIRMED** |
| Terminal velocity (Beard) | Fortran uses `rho_air(k)`, Python overwrites with `1.0` | **BUG CONFIRMED** |
| Same-weight collision | Fortran uses `FLOOR(A/2)`, Python uses `A*0.5` | Minor difference |

---

## Assessment Team

| Role | Key Findings |
|------|-------------|
| **Science Auditor** | Surface tension formula wrong at warm temps; collision module has fundamental bugs (V_parcel, E_H80 units, kappa); no ventilation/radiation (acceptable for education) |
| **Education Expert** | Score: 5/10 as teaching tool. Zero explanatory text, no guided experiments, no Kohler curve visualization. Strong physics engine needs pedagogical wrapper |
| **UI Developer** | Plotly memory leak is critical. No input validation beyond widget bounds. Jet colormap, deepcopy-in-loop, no progress bar. FigureWidget switch is biggest quick win |
| **Code Reviewer** | 6 CRITICAL, 6 HIGH, 8 MEDIUM, 8 LOW issues. No tests. 38-value function returns. Wildcard imports everywhere. All collision results quantitatively wrong |

---

## Files Changed Summary (Estimated)

| File | Changes Needed | Priority |
|------|---------------|----------|
| `condensation.py` | Fix surface tension, remove dead V_parcel, remove unused numba import | P1 |
| `collision.py` | Fix V_parcel/rho_parcel, fix E_H80 units, fix kappa corruption, remove hardcoded constants | P1 |
| `parameters.py` | Fix N_l tuple, fix rho_a units | P1 |
| `analysis.py` | Fix division by zero | P1 |
| `parcel.py` | Fix temperature sign logic | P1 |
| `aero_init.py` | Use esatw() consistently, add max_iter to r_equi, fix kappa for Weighting_factor | P1-P2 |
| `PyLCM_edu.ipynb` | Add markdown cells, experiments, Kohler viz, presets | P2 |
| `animation.py` | FigureWidget refactor | P2-P3 |
| `widget.py` | Bounded nt, sigma/mu bounds, presets, stop button | P2-P3 |
| `print_plot.py` | Viridis colormap, fix typo, fix deepcopy, fix deprecated API | P3 |
| `timestep_routine.py` | Progress bar, move time_array, remove unused imports | P3 |
