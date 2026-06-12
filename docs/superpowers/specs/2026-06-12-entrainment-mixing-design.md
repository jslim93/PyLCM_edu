# Entrainment / Mixing (Warm-Cloud, Parameterized + LEM-ready) — Design Spec

**Date:** 2026-06-12
**Author:** J. Lim (jslim93) with agent assistance
**Base branch:** `feature/performance-ensemble` (entrainment branches off it directly)
**Status:** Approved design — pending user spec review → implementation plan

---

## 1. Motivation & scope

PyLCM's entrainment is currently a crude homogeneous bulk dilution (`basic_entrainment`),
hard-gated as experimental. This cycle makes entrainment a real, teachable feature: the
**homogeneous-vs-inhomogeneous mixing** contrast, parameterized by the **Inhomogeneous Mixing
Degree (IHMD)** of Lim & Hoffmann (2023, *JGR Atmos.*, "Between Broadening and Narrowing").

**In scope (this cycle):** warm-cloud (liquid) entrainment mixing, parameterized by a mixing
fraction χ and IHMD, with an interface that a Linear Eddy Model (LEM) backend can later
implement without driver changes.

**Out of scope (deferred):**
- **Phase 4 — ice / mixed-phase** (the `LEM_MIXED_PHASE_CONCEPTS.md` hard problems:
  sedimentation in the mixing domain, ice-particle statistics, Wegener–Bergeron–Findeisen).
- A full LEM implementation (the interface is provided; the backend is a later cycle).

## 2. Physics

**Environmental profiles.** `create_env_profiles` (parcel.py) supplies `T_env(z)`, `q_env(z)`.

**Mixing fraction χ.** For a conservative scalar X, χ = (X − X_env)/(X_cloud − X_env): the
fraction of cloudy air in the mixture (χ=1 pure cloud, χ=0 pure environment). [FACT, vault]
Entrainment progressively lowers χ as dry air is mixed in.

**Entrainment rate λ [m⁻¹].** Per step the parcel entrains a dry-air mass fraction
`δ = λ · w · dt` (= λ·dz). Bulk conservative scalars relax toward the environment:
`q ← q + δ(q_env − q)`, `T ← T + δ(T_env − T)`.

**IHMD mechanism (the teaching centerpiece).** The bulk dilution lowers saturation, so a
liquid-water amount ΔLWC must evaporate to re-equilibrate. IHMD ∈ [0,1] controls **how that
evaporation is distributed across super-droplets**, targeting the confirmed relation:

> **N_c / N_{c,0} = (q_c / q_{c,0})^IHMD**   (IHMD=0 homogeneous, IHMD=1 inhomogeneous)

Endpoint behavior (both [FACT, vault] / author-confirmed):
- **IHMD = 0 (homogeneous):** distribute ΔLWC across *all* droplets (each shrinks); N conserved;
  DSD broadens toward small sizes.
- **IHMD = 1 (inhomogeneous):** evaporate *whole* super-droplets (reduce multiplicity/remove,
  returning water to vapor and the core to aerosol) until ΔLWC is met; survivors keep size; N drops.
- **0 < IHMD < 1:** split ΔLWC between the two mechanisms so the N–q power law above holds.

**Water budget.** Evaporated liquid returns to vapor; total water (vapor + liquid + aerosol
solute) is conserved by construction at every mixing step.

## 3. Components & interfaces (LEM-ready)

- `PyLCM/mixing.py` — new module. A `MixingModel` protocol:
  `apply(particles_list, T, q, P, z, dt, w, rng) -> (particles_list, T, q)`.
- `ParameterizedMixing(lambda_ent, ihmd, env_profiles)` — implements §2 (this cycle).
- `LEMMixing(...)` — **interface stub only** this cycle (raises `NotImplementedError` with a
  pointer to Phase-3b); proves the seam is real. A future cycle fills it (1D triplet-map domain,
  per-droplet S′; Hoffmann & Feingold 2019 style).
- The timestep driver only ever calls `MixingModel.apply(...)`, so swapping in LEM later needs
  no driver edits.

The existing `entrainment.py` `basic_entrainment` is superseded by `ParameterizedMixing`
(homogeneous limit, IHMD=0) and retired or left gated; no other module imports its internals.

## 4. Integration & controls

- **Per-timestep order: mixing first, then condensation.** ascent → `MixingModel.apply(...)` runs
  as one self-contained step (dilute bulk `T`,`q` toward the environment; compute the ΔLWC that
  must evaporate; distribute it as *homogeneous* shrink of all droplets and *inhomogeneous* whole-
  droplet removal so that `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD`; return evaporated water to vapor) →
  the existing condensation step then performs normal diffusional growth/evaporation on the mixed
  state → collision. Mixing owns the entrainment-evaporation; condensation is unchanged.
- **Two widget sliders** (teaching UX): `entrainment_rate` (λ) and `mixing_degree` (**IHMD**,
  0=homogeneous → 1=inhomogeneous). χ is diagnosed from the dilution, not set directly.
- Default `entrainment_enabled = False` so all existing runs/notebooks are byte-for-byte unchanged.
- Inhomogeneous removal uses `np.random` (seedable) — consistent with the rest of the model.

## 5. Testing strategy

Physics invariants (survive Monte-Carlo variance):
- **Water budget:** vapor + liquid + solute conserved across a mixing step to tight tolerance.
- **Homogeneous (IHMD=0):** droplet number N_c conserved; mean volume radius r_v decreases.
- **Inhomogeneous (IHMD=1):** N_c decreases; surviving-droplet mean size ~unchanged (within tol).
- **N–q power law:** after a mixing step, `N_c/N_{c,0} ≈ (q_c/q_{c,0})^IHMD` across IHMD ∈ {0,0.5,1}.
- **Determinism:** fixed `np.random.seed` ⇒ identical result.
- **No-op:** `entrainment_enabled=False` reproduces the pre-feature golden run bit-for-bit.

## 6. Validation deliverable

`validation/entrainment_mixing.ipynb`: identical ascent run at fixed λ, sweeping IHMD ∈ {0, 0.5, 1};
plot the resulting DSDs and the N–r_v mixing diagram — reproducing the homogeneous-broadening vs
inhomogeneous-number-depletion contrast of Lim & Hoffmann (2023).

## 7. Branch / workspace

New git worktree `../PyLCM_entrain` on branch `feature/entrainment-mixing`, off
`feature/performance-ensemble` (so validation runs on the ~3× faster engine). Never touch the
live `/Users/dr.cloud/PyLCM_edu` session. No "Claude"/Co-Authored-By commit lines.

## 8. Acceptance criteria

- `ParameterizedMixing` passes all §5 invariants in CI.
- The `MixingModel` seam exists and `LEMMixing` stub raises with a clear Phase-3b pointer.
- Validation notebook renders the IHMD sweep (DSD + N–r_v diagram).
- `entrainment_enabled=False` leaves existing behavior unchanged (no-op golden test passes).
- Two teaching sliders (λ, IHMD) wired through widget → driver → mixing module.

## 9. Roadmap context

- **Phase 3 (this spec):** warm-cloud parameterized IHMD mixing + LEM interface.
- **Phase 3b (later):** LEM backend behind the same interface (1D triplet map, per-droplet S′).
- **Phase 4 (later):** ice / mixed-phase mixing (`LEM_MIXED_PHASE_CONCEPTS.md`).
